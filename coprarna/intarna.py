"""IntaRNAsTar orchestration, compressed interactions and artifact filtering."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import json
import hashlib
import re
import tempfile
from Bio import SeqIO
from . import genome, tables, phylogeny, energy, parameters
from .io import accession, open_text, read_text, write_text
from .runtime import options, command


def numeric_sort(lines, index=0, separator=";"):
    """Reproduce numeric table sorting with deterministic text ties and header first."""

    def key(line):
        fields = line.rstrip("\n").split(separator)
        try:
            return (1, float(fields[index]), line)
        except (ValueError, IndexError):
            return (0, 0, line)

    return sorted(lines, key=key)


def perfect_helix(row, windows=None, fraction=0.9):
    """A gap-free helix spanning >=90% of either RNA's effective window."""
    structure = row["hybridDP"].split("&")
    if len(structure) != 2:
        return False
    left, right = structure
    if not left or set(left) != {"("} or set(right) != {")"} or len(left) != len(right):
        return False
    windows = windows or {}
    for key, sequence in [("tAccW", "seq1"), ("qAccW", "seq2")]:
        length = len(row[sequence])
        window = int(windows.get(key, 0))
        effective = min(length, window) if window > 0 else length
        if effective and len(left) >= fraction * effective:
            return True
    return False


def filter_perfect_helices(path, windows=None):
    """Atomically filter plain/gzip native CSV and return an auditable removal list.

    Read one interaction at a time. The filtered table replaces its source before
    optimal/suboptimal splitting or GEV fitting, so excluded energies cannot affect
    the estimated null distribution. Keep identifiers and coordinates in the audit."""
    path = Path(path)
    removed = []
    kept = 0
    with tempfile.NamedTemporaryFile(
        prefix=".filtered-", suffix=path.suffix, dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
    try:
        with open_text(path) as source, open_text(temporary, "wt") as target:
            reader = csv.DictReader(source, delimiter=";")
            writer = csv.DictWriter(
                target, fieldnames=reader.fieldnames, delimiter=";", lineterminator="\n"
            )
            writer.writeheader()
            for row in reader:
                if perfect_helix(row, windows):
                    removed.append(
                        {
                            key: row[key]
                            for key in (
                                "id1",
                                "id2",
                                "start1",
                                "end1",
                                "start2",
                                "end2",
                            )
                        }
                    )
                else:
                    writer.writerow(row)
                    kept += 1
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return {
        "removed_perfect_helices": len(removed),
        "retained_interactions": kept,
        "removed_interactions": removed,
    }


def prediction_signature(directory, target, query, *, native_version=None):
    """Hash scientific inputs and settings to prevent stale prediction reuse.

    Include tool version, personality, filter policy and file contents. Exclude
    absolute paths so equivalent jobs and inline/file settings have the same key.
    Native callers supply the detected release instead of trusting stale options."""
    digest = hashlib.sha256(b"CopraRNA2:perfect-helix-filter=.9:v1")
    version = (
        native_version
        if native_version is not None
        else options(directory).get("IntaRNA version", "")
    )
    digest.update(version.encode())
    digest.update(parameters.PERSONALITY.encode())
    files = [Path(directory) / target, Path(directory) / query]
    params = parameters.arguments(directory)
    if params:
        files.append(Path(params[1]))
    for path in files:
        digest.update(b"\0")
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def prepare_predictions(
    srna: Path, up: int, down: int, region: str, groups: list[str], directory: Path
) -> str:
    """Generate target/query FASTAs, predict and cluster compressed interactions.

    Groups contain comma-separated replicon paths, with the primary chromosome
    first. Each organism runs independently using the same parameter file. Only
    matching cached predictions are reused; filtering precedes all statistical work."""
    native_version = parameters.check_version(directory)
    opt = options(directory)
    for group in groups:
        main = accession(group.split(",")[0])
        target = directory / f"{main}_upfromstartpos_{up}_down_{down}.fa"
        with target.open("w") as stream:
            for filename in group.split(","):
                stream.write(
                    genome.parse_regions(directory / filename, up, down, region)
                )
    queries = {}
    for record in SeqIO.parse(srna, "fasta"):
        queries[record.id.removeprefix("ncRNA_")] = record.id + ".fa"
        (directory / (record.id + ".fa")).write_text(
            f">{record.id}\n" + re.sub(r"\s|-", "", str(record.seq)).lower()
        )
    windows = parameters.values(directory)

    def predict(group):
        main = accession(group.split(",")[0])
        target = f"{main}_upfromstartpos_{up}_down_{down}.fa"
        output = target + ".intarna.csv.gz"
        signature = prediction_signature(
            directory, target, queries[main], native_version=native_version
        )
        report_path = directory / (target + ".filter.json")
        try:
            cached = json.loads(report_path.read_text())
        except (OSError, ValueError):
            cached = {}
        if (
            not (directory / output).exists()
            or cached.get("native_signature") != signature
        ):
            command(
                [
                    parameters.PERSONALITY,
                    *parameters.arguments(directory),
                    *([] if "outNumber" in windows else ["--outNumber", "2"]),
                    "--target",
                    target,
                    "--query",
                    queries[main],
                    "--outCsvCols",
                    tables.INTARNA_HEADER.replace(";", ","),
                    "--outMode=C",
                    "--out",
                    output,
                ],
                directory,
            )
            if not (directory / output).is_file():
                raise RuntimeError(f"IntaRNAsTar did not produce {output}")
            report = filter_perfect_helices(directory / output, windows)
            report.update(
                native_signature=signature,
                personality=parameters.PERSONALITY,
                perfect_helix_fraction=0.9,
            )
            report_path.write_text(json.dumps(report, indent=2) + "\n")
        sorted_path = directory / (target + ".intarna.sorted.csv.gz")
        write_text(
            sorted_path,
            "".join(
                numeric_sort(
                    read_text(directory / output).splitlines(keepends=True), 14
                )
            ),
        )

    with ThreadPoolExecutor(max_workers=max(1, int(opt["core count"]))) as pool:
        list(pool.map(predict, groups))
    tables.disentangle(directory)
    energy.add_energy_pvalues(directory)
    (directory / "opt_tags.clustered").write_text(tables.cluster_predictions(directory))
    command(
        ["mafft", "--localpair", "--quiet", "16s_sequences.fa"],
        directory,
        output="16s_sequences.aln",
    )
    command(
        [
            "distmat",
            "-sequence",
            "16s_sequences.aln",
            "-nucmethod",
            "1",
            "-outfile",
            "distmat.out",
            "-auto",
        ],
        directory,
    )
    (directory / "compatible.distmat").write_text(
        phylogeny.transform_distance(directory / "distmat.out")
    )
    return ""
