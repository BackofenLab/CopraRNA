"""Every migrated fixture must match independently captured legacy output."""

import json
from pathlib import Path
import shutil

import pytest

from coprarna.helpers import HELPERS, run_helper
from tests.support.legacy import FIXTURES, canonicalize

MIGRATED = {helper.legacy for helper in HELPERS.values()}
CASES = [
    p.parent
    for p in sorted(FIXTURES.glob("*/case.json"))
    if Path(json.loads(p.read_text())["script"]).name in MIGRATED
]


@pytest.mark.parametrize("case", CASES, ids=[p.name for p in CASES])
def test_python_matches_golden_bytes(case, tmp_path):
    spec = json.loads((case / "case.json").read_text())
    shutil.copytree(case / "input", tmp_path, dirs_exist_ok=True)
    stdout = run_helper(
        Path(spec["script"]).name, spec["args"], tmp_path, seed=spec.get("rng_seed")
    )
    actual = {"stdout": stdout.encode()}
    actual.update({name: (tmp_path / name).read_bytes() for name in spec["outputs"]})
    for name, value in actual.items():
        expected = (case / "expected" / name).read_bytes()
        assert canonicalize(name, value, spec) == canonicalize(name, expected, spec)
