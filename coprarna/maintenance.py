"""Organism mapping maintenance, including the original pseudo-KEGG identifiers."""

from datetime import datetime
from pathlib import Path
import random
import re
from urllib.request import urlopen


class PerlRandom:
    def __init__(self, seed=None):
        if seed is None:
            seed = random.SystemRandom().getrandbits(32)
        self.state = ((int(seed) << 16) + 0x330E) & ((1 << 48) - 1)

    def random(self):
        self.state = (0x5DEECE66D * self.state + 11) & ((1 << 48) - 1)
        return self.state / (1 << 48)


def write_availability_header(directory: Path, header: str, now=None) -> str:
    now = now or datetime.now().astimezone()
    date = now.strftime("%a %b %e %H:%M:%S %Z %Y")
    text = (
        "## last updated: "
        + date
        + "\n"
        + header
        + (directory / "CopraRNA_available_organisms.tmp").read_text()
    )
    (directory / "CopraRNA_available_organisms.txt").write_text(text)
    (directory / "CopraRNA_available_organisms.tmp").unlink()
    return ""


def build_organism_mapping(
    directory: Path, report: Path | None = None, seed=None, now=None
) -> str:
    if report is None:
        report = directory / "prokaryotes.txt"
        with urlopen(
            "ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
        ) as response:
            report.write_bytes(response.read())
    rng = PerlRandom(seed)
    used = set()
    printed = set()
    mapping = []
    available = []
    for line in report.read_text().splitlines():
        if not re.search("N[CZ]_.+", line):
            continue
        code = ""
        while not code or code in used:
            code = "".join(chr(97 + int(rng.random() * 26)) for _ in range(4))
        used.add(code)
        fields = line.split("\t")
        ids = []
        for entry in fields[8].split(";"):
            match = re.search("(N[CZ]_.+)", entry)
            if match:
                ids.append(re.sub(r"\.\d+", "", match[1]).split("/")[0])
        name = re.sub(r"\s+", "_", fields[0].strip())
        help_line = " ".join(ids) + "\t" + name + "\n"
        if help_line not in printed:
            available.append(help_line)
            mapping.append(code + "\t" + " ".join(ids) + "\n")
        printed.add(help_line)
    (directory / "kegg2refseqnew.csv").write_text("".join(mapping))
    (directory / "CopraRNA_available_organisms.tmp").write_text("".join(available))
    from importlib.resources import files

    header = files("coprarna").joinpath("data/organism_head.txt").read_text()
    return write_availability_header(directory, header, now)
