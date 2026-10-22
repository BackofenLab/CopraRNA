import gzip
import importlib
from importlib.resources import files
import json
from pathlib import Path

from coprarna.helpers import HELPERS
from tests.support.legacy import ROOT


def test_mapping_covers_every_original_script_and_discloses_pending_work():
    original = json.loads((ROOT / "tests/fixtures/legacy-sha256.json").read_text())
    scripts = {
        name
        for name in original
        if Path(name).suffix in {".pl", ".R", ".py", ".js", ".sh"}
    }
    inventory = json.loads(files("coprarna").joinpath("migration.json").read_text())
    assert {entry["legacy"] for entry in inventory} == scripts
    assert len(inventory) == len(scripts)
    ported = {
        Path(entry["legacy"]).name for entry in inventory if entry["status"] == "ported"
    }
    assert ported == {helper.legacy for helper in HELPERS.values()}
    for entry in inventory:
        mapping = (ROOT / (entry["legacy"] + ".md")).read_text()
        assert entry["status"] in mapping
        if entry["status"] == "ported":
            assert entry["fixtures"]
            module, function = entry["target"].rsplit(".", 1)
            assert callable(getattr(importlib.import_module(module), function))
        else:
            assert entry["status"] == "retired"
            assert "Retired" in mapping
            assert "target" not in entry


def test_packaged_organism_tables_are_unchanged():
    for name in ("kegg2refseqnew.csv", "CopraRNA_available_organisms.txt"):
        bundled = files("coprarna").joinpath("data/" + name + ".gz").read_bytes()
        assert gzip.decompress(bundled) == (ROOT / "coprarna_aux" / name).read_bytes()
