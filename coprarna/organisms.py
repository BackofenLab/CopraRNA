"""Resolve RefSeq assemblies on demand, with reusable, versioned metadata."""

from datetime import datetime, timezone
import json
from pathlib import Path
import re
import tempfile
from urllib.parse import urlencode

from . import ncbi

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/"
SEQUENCE = r"(?:NC_\d+|NZ_[A-Z]*\d+)"
ASSEMBLY = r"GCF_\d+"


def validate_accession(identifier):
    """Return an NC_/NZ_ sequence or GCF_ assembly ID, optionally versioned.

    Reject malformed identifiers before network access or cache-path construction.
    Only RefSeq genomic identifiers are accepted; this performs no file writes.
    """
    if not isinstance(identifier, str) or not re.fullmatch(
        rf"(?:{SEQUENCE}|{ASSEMBLY})(?:\.[1-9]\d*)?", identifier
    ):
        raise ValueError(
            f"Invalid RefSeq accession: {identifier!r}; use NC_, NZ_ or GCF_"
        )
    return identifier


def _json(url, opener):
    with opener(url) as response:
        result = json.load(response)
    if not isinstance(result, dict) or result.get("error") or result.get("errors"):
        raise ValueError(f"NCBI metadata request failed: {result}")
    return result


def _reports(endpoint, opener):
    """Read every report page; never silently discard chromosomes or plasmids."""
    result, tokens = [], set()
    token = ""
    while True:
        url = (
            DATASETS
            + endpoint
            + "?"
            + urlencode({"page_size": 1000, **({"page_token": token} if token else {})})
        )
        page = _json(url, opener)
        result.extend(page.get("reports", []))
        token = page.get("next_page_token", "")
        if not token:
            if "total_count" in page and len(result) != int(page["total_count"]):
                raise ValueError("Incomplete NCBI assembly report")
            return result
        if token in tokens:
            raise ValueError("Repeated page token in NCBI assembly report")
        tokens.add(token)


def _fetch(identifier, opener):
    if identifier.startswith("GCF_"):
        candidates = [identifier]
        versioned = None
    else:
        summary = _json(
            ncbi.BASE
            + "esummary.fcgi?"
            + urlencode(
                {
                    "db": "nuccore",
                    "id": identifier,
                    "retmode": "json",
                    "tool": "coprarna",
                }
            ),
            opener,
        ).get("result", {})
        records = [summary[uid] for uid in summary.get("uids", [])]
        if len(records) != 1 or not records[0].get("accessionversion"):
            raise ValueError(f"No unique RefSeq sequence found for {identifier}")
        versioned = records[0]["accessionversion"]
        # Datasets sequence_assemblies requires the version even for current IDs.
        if versioned.split(".")[0] != identifier.split(".")[0] or (
            "." in identifier and versioned != identifier
        ):
            raise ValueError(f"NCBI returned {versioned} for {identifier}")
        candidates = _json(
            DATASETS + f"sequence_accession/{versioned}/sequence_assemblies", opener
        ).get("accessions", [])
        candidates = [a for a in candidates if re.fullmatch(ASSEMBLY + r"\.\d+", a)]
    if not candidates:
        raise ValueError(f"No RefSeq assembly found for {identifier}")
    reports = _reports("accession/" + ",".join(candidates) + "/dataset_report", opener)
    current = {
        r["accession"]: r
        for r in reports
        if (
            re.fullmatch(ASSEMBLY + r"\.\d+", r.get("accession", ""))
            and r.get("assembly_info", {}).get("assembly_status") == "current"
            and r["assembly_info"].get("assembly_level") == "Complete Genome"
        )
    }
    if len(current) != 1:
        raise ValueError(
            f"Expected one current complete RefSeq assembly for {identifier}; found "
            f"{', '.join(sorted(current)) or 'none'}. "
            "For ambiguous sequences, use a specific GCF_ accession in the FASTA header."
        )
    assembly, report = next(iter(current.items()))
    sequences = _reports(f"accession/{assembly}/sequence_reports", opener)
    if not sequences or any(
        s.get("role") != "assembled-molecule"
        or not re.fullmatch(SEQUENCE + r"\.\d+", s.get("refseq_accession", ""))
        for s in sequences
    ):
        raise ValueError(f"{assembly} does not contain only complete RefSeq replicons")
    # One chromosome labels each organism throughout the scientific workflow.
    # Longest chromosome first, then other chromosomes and plasmids, deterministically.
    sequences.sort(
        key=lambda s: (
            s.get("assigned_molecule_location_type") != "Chromosome",
            -int(s.get("length", 0)),
            s["refseq_accession"],
        )
    )
    replicons = [s["refseq_accession"] for s in sequences]
    if versioned and versioned not in replicons:
        raise ValueError(
            f"{versioned} is not a replicon of current assembly {assembly}"
        )
    return {
        "schema_version": 1,
        "query": identifier,
        "assembly": assembly,
        "organism": report["organism"]["organism_name"],
        "replicons": replicons,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }


def validate_metadata(metadata, identifier):
    """Validate a cached resolution and return it without changing files.

    Require complete versioned replicon identities and membership of the requested
    sequence/assembly. Malformed or mismatched metadata raises ValueError.
    """
    validate_accession(identifier)
    if not isinstance(metadata, dict):
        raise ValueError(f"Invalid genome metadata for {identifier}")
    replicons = metadata.get("replicons")
    if (
        metadata.get("schema_version") != 1
        or metadata.get("query") != identifier
        or not re.fullmatch(ASSEMBLY + r"\.\d+", str(metadata.get("assembly", "")))
        or not isinstance(metadata.get("organism"), str)
        or not metadata["organism"]
        or not isinstance(replicons, list)
        or not replicons
        or any(
            not isinstance(r, str) or not re.fullmatch(SEQUENCE + r"\.\d+", r)
            for r in replicons
        )
    ):
        raise ValueError(f"Invalid genome metadata for {identifier}")
    if len({r.split(".")[0] for r in replicons}) != len(replicons):
        raise ValueError(f"Duplicate replicons in genome metadata for {identifier}")
    members = [metadata["assembly"]] if identifier.startswith("GCF_") else replicons
    if identifier not in (
        members if "." in identifier else [r.split(".")[0] for r in members]
    ):
        raise ValueError(f"Genome metadata does not contain {identifier}")
    return metadata


def resolve(identifier, database, *, refresh=False, opener=None):
    """Resolve one RefSeq ID to its assembly, organism name and ordered replicons.

    Cache JSON under database/metadata/<query>.json using an accession lock and
    atomic replacement. Cached resolutions permit offline use and stay fixed until
    refresh=True. Failed refreshes preserve the previous entry. Genome sequences
    are downloaded separately by ncbi.genome_path; opener is injectable for tests.
    """
    validate_accession(identifier)
    directory = Path(database).expanduser().resolve() / "metadata"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (identifier + ".json")
    with ncbi._lock(directory / (identifier + ".lock")):
        if path.exists() and not refresh:
            return validate_metadata(json.loads(path.read_text()), identifier)
        metadata = validate_metadata(
            _fetch(identifier, opener or ncbi.open_url), identifier
        )
        with tempfile.NamedTemporaryFile(
            mode="w", dir=directory, delete=False
        ) as stream:
            temporary = Path(stream.name)
            try:
                json.dump(metadata, stream, indent=2)
                stream.write("\n")
                stream.close()
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
    return metadata


def resolve_many(identifiers, database, *, refresh=False):
    """Resolve an ordered job's IDs and reject duplicate/overlapping assemblies.

    Populate the shared metadata cache but do not download genomes or write job
    outputs. Different chromosomes/plasmids of one assembly cannot count as
    independent organisms in comparative statistics.
    """
    result, assemblies, replicons = [], set(), set()
    for identifier in identifiers:
        metadata = resolve(identifier, database, refresh=refresh)
        # A refreshed alias may point to a newer version than another cached ID.
        # Sequence/assembly revisions still represent the same organism.
        assembly = metadata["assembly"].split(".")[0]
        members = {r.split(".")[0] for r in metadata["replicons"]}
        if assembly in assemblies or replicons.intersection(members):
            raise ValueError(f"Duplicate organism/assembly in input: {identifier}")
        assemblies.add(assembly)
        replicons.update(members)
        result.append(metadata)
    return result
