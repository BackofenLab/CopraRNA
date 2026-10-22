"""Every migrated fixture must match independently captured legacy output."""

import json
import shutil

import pytest

from coprarna.helpers import run_helper
from tests.support.fixtures import FIXTURES, canonicalize

# Retain independently captured bytes. Cases whose contract changed have explicit
# assertions below or dedicated bugfix tests, never regenerated golden files.
CASES = [p.parent for p in sorted(FIXTURES.glob("*/case.json"))]


@pytest.mark.parametrize("case", CASES, ids=[p.name for p in CASES])
def test_python_matches_golden_bytes(case, tmp_path):
    spec = json.loads((case / "case.json").read_text())
    shutil.copytree(case / "input", tmp_path, dirs_exist_ok=True)
    if spec["helper"] == "extract-proteins":
        from coprarna.genome import extract_proteins

        # Code assignment is now explicit; retain the captured translated sequences.
        stdout = extract_proteins(tmp_path / spec["args"][0], organism="esgb")
    else:
        stdout = run_helper(spec["helper"], spec["args"], tmp_path)
    actual = {"stdout": stdout.encode()}
    actual.update({name: (tmp_path / name).read_bytes() for name in spec["outputs"]})
    for name, value in actual.items():
        if spec["helper"] == "disentangle" and "subopt" in name:
            # The corrected second-best contract is asserted in test_bugfixes.
            continue
        expected = (case / "expected" / name).read_bytes()
        assert canonicalize(name, value, spec) == canonicalize(name, expected, spec)
