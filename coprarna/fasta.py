"""FASTA formatting and DomClust sequence metadata (legacy byte formats)."""

from pathlib import Path
import re

from Bio import SeqIO


def format_fasta(path: Path) -> str:
    """Return normalized FASTA using record IDs and unwrapped sequences.

    Discard descriptions and sequence whitespace while preserving record order
    and base case. Read the source without changing or moving it.
    """
    output = []
    for record in SeqIO.parse(path, "fasta"):
        sequence = re.sub(r"\s", "", str(record.seq))
        output.append(f">{record.id}\n{sequence}\n")
    return "".join(output)


def fasta_to_genefile(path: Path) -> str:
    """Write DomClust .gene and .tit metadata beside a protein FASTA.

    Interpret identifiers as organism:gene. Preserve input order and emit protein
    lengths plus descriptions in the retained text format; return empty stdout.
    """
    genes, titles = [], []
    name, title, sequence = "", "", ""

    def emit():
        if sequence:
            fields = name.split(":")
            species = fields[0]
            gene = fields[1] if len(fields) > 1 else ""
            genes.append(f"{species} {gene} {len(sequence)}\n")
            titles.append(f"{species}:{gene}\t{title}\n")

    for line in path.read_text().splitlines():
        match = re.match(r"^>\s*(\S.*)", line)
        if match:
            emit()
            parts = re.split(r"\s+", match[1], maxsplit=1)
            name, title = parts[0], parts[1] if len(parts) > 1 else ""
            sequence = ""
        else:
            sequence += re.sub(r"\s", "", line)
    emit()
    Path(str(path) + ".gene").write_text("".join(genes))
    Path(str(path) + ".tit").write_text("".join(titles))
    return ""
