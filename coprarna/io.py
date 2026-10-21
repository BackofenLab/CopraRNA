"""Streaming access to compressed and plain scientific data files."""

from contextlib import contextmanager
import gzip
from pathlib import Path
from Bio import SeqIO


@contextmanager
def open_text(path, mode="rt"):
    """Open UTF-8 text directly from plain or gzip storage; preserve CSV newlines."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, mode, encoding="utf-8", newline="") as stream:
        yield stream


def read_text(path):
    """Read a small complete text table through the shared compression adapter."""
    with open_text(path) as stream:
        return stream.read()


def write_text(path, text):
    """Write a text table as gzip when its destination has a .gz suffix."""
    with open_text(path, "wt") as stream:
        stream.write(text)


def paths(directory, pattern):
    """Prefer compressed data when both representations of a table exist."""
    files = {p.name: p for p in Path(directory).glob(pattern)}
    for path in Path(directory).glob(pattern + ".gz"):
        files[path.name[:-3]] = path
    return [files[name] for name in sorted(files)]


def accession(path):
    """Extract the RefSeq label from plain/compressed GenBank file names."""
    name = Path(path).name.removesuffix(".gz")
    for suffix in (".gbff", ".gbk", ".gb"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def genbank_records(path):
    """Yield GenBank records from the compressed stream without a plain disk copy."""
    with open_text(path) as stream:
        yield from SeqIO.parse(stream, "genbank")


def genome_paths(directory):
    """Read the job genome manifest, or discover local genomes for standalone helpers."""
    manifest = Path(directory) / "genome_paths.txt"
    if manifest.exists():
        return [Path(line) for line in manifest.read_text().splitlines() if line]
    return sorted(
        {
            p
            for suffix in ("*.gb", "*.gb.gz", "*.gbk.gz", "*.gbff.gz")
            for p in Path(directory).glob(suffix)
        }
    )
