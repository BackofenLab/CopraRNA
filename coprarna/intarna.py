"""Native orchestration of whole-genome IntaRNA predictions and phylogeny input."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
from Bio import SeqIO
from . import genome, tables, phylogeny, energy, statistics
from .runtime import options, command


def numeric_sort(lines, index=0, separator=";"):
    def key(line):
        fields = line.rstrip("\n").split(separator)
        try:
            number = float(fields[index])
            return (1, number, line)
        except (ValueError, IndexError):
            return (0, 0, line)

    return sorted(lines, key=key)


def prepare_predictions(
    srna: Path, up: int, down: int, region: str, groups: list[str], directory: Path
) -> str:
    opt = options(directory)
    for group in groups:
        main = group.split(",")[0][:-3]
        output = f"{main}_upfromstartpos_{up}_down_{down}.fa"
        for i, filename in enumerate(group.split(","), 1):
            (directory / f"{output}{i}").write_text(
                genome.parse_regions(directory / filename, up, down, region)
            )
    filenames = sorted(path.name for path in directory.iterdir())
    for group in groups:
        main = group.split(",")[0][:-3]
        output = f"{main}_upfromstartpos_{up}_down_{down}.fa"
        with (directory / output).open("a") as stream:
            for name in filenames:
                if re.search(re.escape(main) + r"\S+\.fa\d+", name):
                    stream.write((directory / name).read_text())
    queries = {}
    for record in SeqIO.parse(srna, "fasta"):
        name = record.id
        queries[name.removeprefix("ncRNA_")] = name + ".fa"
        (directory / (name + ".fa")).write_text(
            f">{name}\n" + re.sub(r"\s|-", "", str(record.seq)).lower()
        )

    def predict(group):
        main = group.split(",")[0][:-3]
        target = f"{main}_upfromstartpos_{up}_down_{down}.fa"
        output = target + ".intarna.csv"
        if not (directory / output).exists():
            command(
                [
                    "IntaRNA",
                    "--tAccW",
                    opt["win size"],
                    "--tAccL",
                    opt["max bp dist"],
                    "--outNumber",
                    "2",
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

    with ThreadPoolExecutor(max_workers=max(1, int(opt["core count"]))) as pool:
        list(pool.map(predict, groups))
    for path in sorted(directory.glob("*intarna.csv")):
        (directory / (path.name[:-3] + "sorted.csv")).write_text(
            "".join(numeric_sort(path.read_text().splitlines(keepends=True), 14))
        )
    tables.disentangle(directory)
    energy.add_energy_pvalues(directory)
    (directory / "opt_tags.clustered").write_text(tables.cluster_predictions(directory))
    if int(opt["CopraRNA1"]):
        statistics.filter_clusters(directory)
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
