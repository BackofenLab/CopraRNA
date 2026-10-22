"""Refresh bundled reference genomes from NCBI and record download provenance."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tempfile

from coprarna import ncbi, organisms
from coprarna.io import genbank_records


def refresh(directory, accessions):
    """Download current versions into staging, validate records, then publish.

    Replace only the reference GenBank archives and metadata listed in this run.
    Write manifest.json with exact IDs, record dates, lengths and gzip SHA-256s.
    Network/parse failures before publication leave the existing archive intact.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".refresh-", dir=directory) as staging:
        staging = Path(staging)
        entries = []
        for query in accessions:
            metadata = organisms.resolve(query, staging)
            version = next(r for r in metadata["replicons"] if r.split(".")[0] == query)
            path = ncbi.genome_path(version, staging, staging)
            records = list(genbank_records(path))
            if len(records) != 1 or records[0].id != version:
                raise ValueError(
                    f"Downloaded GenBank identity does not match {version}"
                )
            record = records[0]
            entry = {
                "query": query,
                "accession": record.id,
                "file": path.name,
                "assembly": metadata["assembly"],
                "organism": record.annotations["organism"],
                "record_date": record.annotations.get("date"),
                "length": len(record),
                "cds_count": sum(f.type == "CDS" for f in record.features),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "source": ncbi.BASE
                + f"efetch.fcgi?db=nuccore&id={version}&rettype=gbwithparts&retmode=text",
            }
            entries.append(entry)
            print(f"Validated {version}: {len(record):,} bp", flush=True)
        manifest = {
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "genomes": entries,
        }
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (directory / "metadata").mkdir(exist_ok=True)
        for entry in entries:
            (staging / entry["file"]).replace(directory / entry["file"])
            name = entry["query"] + ".json"
            (staging / "metadata" / name).replace(directory / "metadata" / name)
            for old in directory.glob(entry["query"] + ".gb.bz2"):
                old.unlink()
            for old in directory.glob(entry["query"] + ".*.gb.gz"):
                if old.name != entry["file"]:
                    old.unlink()
        (staging / "manifest.json").replace(directory / "manifest.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("core_genome_archive"))
    parser.add_argument("accessions", nargs="*")
    args = parser.parse_args()
    accessions = args.accessions or [
        g["query"]
        for g in json.loads((args.directory / "manifest.json").read_text())["genomes"]
    ]
    refresh(args.directory, accessions)
