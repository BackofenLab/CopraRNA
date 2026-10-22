"""Read the unchanged legacy data shipped with the Python wheel."""

import gzip
from importlib.resources import files


def organism_mapping() -> str:
    """Return the packaged organism-code/replicon table as UTF-8 text.

    Read the unchanged gzip resource through importlib so installed wheels behave
    the same as a source checkout; no working-directory file is needed.
    """
    return gzip.decompress(
        files("coprarna").joinpath("data/kegg2refseqnew.csv.gz").read_bytes()
    ).decode()


def archive_readme() -> str:
    """Return the current report-archive README from the installed package resources."""
    return files("coprarna").joinpath("data/archive_README.txt").read_text()
