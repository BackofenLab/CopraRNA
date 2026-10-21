"""Read the unchanged legacy data shipped with the Python wheel."""

import gzip
from importlib.resources import files


def organism_mapping() -> str:
    return gzip.decompress(
        files("coprarna").joinpath("data/kegg2refseqnew.csv.gz").read_bytes()
    ).decode()


def archive_readme() -> str:
    return files("coprarna").joinpath("data/archive_README.txt").read_text()
