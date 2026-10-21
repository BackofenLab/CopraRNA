"""Check committed golden outputs against untouched legacy programs."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from tests.support.legacy import ROOT, FIXTURES, capture, canonicalize, interpreter

CASES = sorted(FIXTURES.glob("*/case.json"))


@pytest.mark.parametrize("case", [p.parent for p in CASES], ids=[p.parent.name for p in CASES])
def test_legacy_matches_captured_baseline(case, tmp_path):
    spec = json.loads((case / "case.json").read_text())
    executable = interpreter(spec["language"])
    if not executable:
        pytest.skip(f"Legacy {spec['language']} runtime unavailable")
    if spec["language"] == "perl" and spec.get("modules"):
        probe = subprocess.run([executable, *[f"-M{m}" for m in spec["modules"]], "-e", "1"],
                               capture_output=True)
        if probe.returncode:
            pytest.skip("Required legacy Perl modules unavailable")
    for name, actual in capture(case, tmp_path).items():
        expected = (case / "expected" / name).read_bytes()
        assert canonicalize(name, actual, spec) == canonicalize(name, expected, spec)


def test_original_sources_are_unchanged():
    manifest = json.loads((FIXTURES / "legacy-sha256.json").read_text())
    for relative, digest in manifest.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest, relative
