"""Parallel MAFFT target alignments with the legacy file and sequence ordering."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
from .runtime import command, options
from .tables import perl_split


def align_targets(path: Path, directory: Path) -> str:
    """Create ordered target FASTAs and run independent MAFFT alignments in parallel."""
    sequences = {}
    for fasta in sorted(directory.glob("*upfromstartpos*down*fa")):
        lines = fasta.read_text().splitlines(keepends=True)
        for i in range(0, len(lines), 2):
            match = re.search(">(.+)", lines[i])
            if match:
                sequences[match[1].lower()] = lines[i + 1]
    target = directory / "target_alignments"
    target.mkdir(exist_ok=True)
    for i, line in enumerate(path.read_text().splitlines()[1:], 1):
        cells = perl_split(",", line)[2:-3]
        with (target / f"{i}.fa").open("w") as stream:
            for cell in cells:
                match = re.search(r"(.+)\(.+", cell)
                if match:
                    tag = match[1].lower()
                    stream.write(f">{tag}\n{sequences.get(tag, '')}")

    def align(path):
        relative = path.relative_to(directory)
        command(
            ["mafft", "--localpair", "--quiet", str(relative)],
            directory,
            output=str(relative).split(".")[0] + ".aln",
        )

    with ThreadPoolExecutor(
        max_workers=max(1, int(options(directory)["core count"]))
    ) as pool:
        list(pool.map(align, sorted(target.glob("*fa"))))
    return ""
