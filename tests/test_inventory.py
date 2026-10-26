"""Release packaging and independently captured fixture integrity."""

import hashlib
import importlib
import json
import tomllib

from coprarna import __version__
from coprarna.helpers import HELPERS
from tests.support.fixtures import ROOT, FIXTURES


def test_release_version_and_readme_agree():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert project["version"] == __version__ == "3.0.0"
    assert project["readme"] == "README.md"
    assert not (ROOT / "README.2.md").exists()


def test_only_python_helpers_and_current_resources_remain():
    for name in (
        "CopraRNA2.pl",
        "coprarna_aux",
        "update_kegg2refseq",
        "environment-legacy.yml",
        "environment-legacy-r.yml",
        "coprarna/migration.json",
        "coprarna/compatibility.py",
        "coprarna/maintenance.py",
    ):
        assert not (ROOT / name).exists(), name
    for path in (ROOT / "coprarna/data").iterdir():
        assert path.name in {"archive_README.txt", "copra_heatmap.html"}
    for helper in HELPERS.values():
        module, function = helper.target.rsplit(".", 1)
        assert callable(getattr(importlib.import_module(module), function))
    assert "organism-mapping" not in HELPERS and "availability-header" not in HELPERS


def test_captured_scientific_expectations_are_unchanged():
    manifest = json.loads((FIXTURES / "captured-sha256.json").read_text())
    for relative, digest in manifest.items():
        assert (
            hashlib.sha256((FIXTURES / relative).read_bytes()).hexdigest() == digest
        ), relative
