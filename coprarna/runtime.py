"""External scientific executables, options and legacy file operations."""

from pathlib import Path
import os
import subprocess


def options(directory):
    """Read named job options, splitting once so values may contain colons."""
    return dict(
        line.rstrip("\n").split(":", 1)
        for line in (Path(directory) / "CopraRNA_option_file.txt")
        .read_text()
        .splitlines()
        if ":" in line
    )


def command(arguments, directory, output=None, quiet=False):
    """Run a native tool in the job directory with an optional streamed output file.

    Prepend COPRARNA_TOOL_PATH without changing the caller's environment. Surface
    exit failures and IntaRNA's stdout error markers, including releases that print
    an error but return zero. Return decoded stdout when not redirected."""
    environment = os.environ.copy()
    if environment.get("COPRARNA_TOOL_PATH"):
        environment["PATH"] = (
            environment["COPRARNA_TOOL_PATH"] + os.pathsep + environment.get("PATH", "")
        )
    stream = (Path(directory) / output).open("wb") if output else subprocess.PIPE
    try:
        result = subprocess.run(
            [str(arg) for arg in arguments],
            cwd=directory,
            env=environment,
            stdout=stream,
            stderr=subprocess.PIPE,
            check=False,
        )
        native_error = (
            (result.stdout or b"").decode(errors="replace") if output is None else ""
        )
        if result.returncode or (
            str(arguments[0]).startswith("IntaRNA") and "# ERROR" in native_error
        ):
            raise RuntimeError(
                f"{arguments[0]} exited {result.returncode}: {(result.stderr.decode(errors='replace') + native_error).strip()}"
            )
        return result.stdout.decode() if result.stdout else ""
    finally:
        if output:
            stream.close()
