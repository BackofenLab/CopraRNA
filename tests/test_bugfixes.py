"""Dedicated regressions for each isolated correction on the phase-6 branch."""

from pathlib import Path
import shutil
from coprarna.compatibility import legacy_behavior
from coprarna import statistics
from tests.support.legacy import FIXTURES


def test_filter_renames_pvalue_column_without_overwriting_normalized_energy(tmp_path):
    shutil.copytree(FIXTURES / "filter-clusters/input", tmp_path, dirs_exist_ok=True)
    with legacy_behavior(False):
        statistics.filter_clusters(tmp_path)
    header = (
        (tmp_path / "opt_tags.clustered_rcsize").read_text().splitlines()[0].split(";")
    )
    assert header[34] == "E_norm"
    assert header[35] == "p-value"
