"""Read report resources from source checkouts and installed wheels."""

from importlib.resources import files


def archive_readme() -> str:
    """Return the current report-archive README from the installed package resources."""
    return files("coprarna").joinpath("data/archive_README.txt").read_text()
