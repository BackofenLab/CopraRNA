"""Execute original sources as a regression oracle, never production code."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures"


def interpreter(language: str) -> str | None:
    variable, executable = {
        "perl": ("COPRARNA_LEGACY_PERL", "perl"),
        "r": ("COPRARNA_LEGACY_RSCRIPT", "Rscript"),
    }[language]
    return os.environ.get(variable) or shutil.which(executable)


def capture(case: Path, directory: Path) -> dict[str, bytes]:
    spec = json.loads((case / "case.json").read_text())
    shutil.copytree(case / "input", directory, dirs_exist_ok=True)
    executable = interpreter(spec["language"])
    if not executable:
        raise FileNotFoundError(f"Legacy {spec['language']} interpreter unavailable")
    env = dict(os.environ, LC_ALL="C", TZ="UTC")
    # This fixes traversal variability only; Perl does not guarantee insertion order.
    env.update(PERL_HASH_SEED="0", PERL_PERTURB_KEYS="0")
    source = ROOT / spec["script"]
    if "rng_seed" in spec:
        wrapper = directory / "oracle-seeded.R"
        wrapper.write_text(
            f"set.seed({int(spec['rng_seed'])});source({json.dumps(str(source))})\n"
        )
        source = wrapper
    result = subprocess.run(
        [executable, str(source), *spec.get("legacy_args", spec["args"])],
        cwd=directory,
        env=env,
        capture_output=True,
        timeout=spec.get("timeout", 60),
    )
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace"))
    outputs = {"stdout": result.stdout}
    outputs.update({name: (directory / name).read_bytes() for name in spec["outputs"]})
    return outputs


def canonicalize(name: str, value: bytes, spec: dict) -> bytes:
    if name in spec.get("unordered", []):
        header, *rows = value.splitlines(keepends=True)
        return header + b"".join(sorted(rows))
    return value
