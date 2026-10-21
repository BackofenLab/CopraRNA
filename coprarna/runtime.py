"""External scientific executables, options and legacy file operations."""

from pathlib import Path
import os
import subprocess


def options(directory):
    return dict(
        line.rstrip("\n").split(":", 1)
        for line in (Path(directory) / "CopraRNA_option_file.txt")
        .read_text()
        .splitlines()
        if ":" in line
    )


def command(arguments, directory, output=None, quiet=False):
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
        if result.returncode:
            raise RuntimeError(
                f"{arguments[0]} exited {result.returncode}: {result.stderr.decode(errors='replace').strip()}"
            )
        return result.stdout.decode() if result.stdout else ""
    finally:
        if output:
            stream.close()
