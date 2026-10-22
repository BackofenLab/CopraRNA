"""IntaRNA 3 personality and parameter-file integration shared by every caller."""

from contextlib import contextmanager
from pathlib import Path
import re
import tempfile
from .runtime import command, options
from .io import ensure_input_is_distinct

PERSONALITY = "IntaRNAsTar"


def check_version(directory):
    """Reject unavailable/unrecognized IntaRNA personalities and releases before 3."""
    output = command([PERSONALITY, "--version"], directory)
    match = re.search(r"IntaRNA\s+(\d+)\.(\d+)(?:\.(\d+))?", output)
    if not match or int(match[1]) < 3:
        raise RuntimeError("CopraRNA2 requires IntaRNA version 3 or newer")
    return match[0]


@contextmanager
def parameter_file(directory, filename=None, text=None):
    """Yield the supplied file or a temporary file containing inline settings.

    Preserve a byte-exact snapshot of a supplied file for reproducibility. Inline
    files live for the whole job, including web exports, and are unlinked on both
    success and failure. No supplied native settings means sTar defaults."""
    if filename is not None and text is not None:
        raise ValueError("Supply an IntaRNA parameter file or inline text, not both")
    directory = Path(directory)
    if filename is not None:
        path = Path(filename)
        path = path if path.is_absolute() else directory / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f"IntaRNA parameter file does not exist: {path}")
        ensure_input_is_distinct(
            filename if Path(filename).is_absolute() else directory / filename,
            directory,
            allow_parameter_snapshot=True,
        )
        snapshot = directory / "IntaRNA_parameters.txt"
        if path != snapshot.resolve():
            snapshot.write_bytes(path.read_bytes())
        yield path
    elif text is not None:
        (directory / "IntaRNA_parameters.txt").write_text(text)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".cfg",
            prefix=".intarna-",
            dir=directory,
            delete=False,
        ) as stream:
            stream.write(text)
            path = Path(stream.name)
        try:
            yield path
        finally:
            path.unlink(missing_ok=True)
    else:
        (directory / "IntaRNA_parameters.txt").write_text("# IntaRNAsTar defaults\n")
        yield None


def arguments(directory):
    """Return native parameter-file arguments, falling back to the saved snapshot.

    The fallback lets subsequent helper calls read settings after an inline
    parameter file has been removed at the end of its prediction context."""
    configured = options(directory).get("intarna parameter file", "")
    if not configured:
        return []
    path = Path(configured)
    if not path.exists():
        path = Path(directory) / "IntaRNA_parameters.txt"
    return ["--parameterFile", str(path)]


def values(directory):
    """Accessibility windows for the post-prediction artificial-helix check."""
    result = {"tAccW": "150", "qAccW": "150"}
    configured = options(directory).get("intarna parameter file", "")
    path = (
        Path(configured) if configured else Path(directory) / "IntaRNA_parameters.txt"
    )
    if not path.is_file():
        path = Path(directory) / "IntaRNA_parameters.txt"
    if path.is_file():
        for line in path.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                key, value = line.split("=", 1)
                result[key.strip()] = value.strip()
    return result
