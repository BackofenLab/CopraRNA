"""Optional real-genome oracle checks using the 17 bundled reference archives."""

import bz2
import os
import subprocess

import pytest

from coprarna.genome import parse_regions
from tests.support.legacy import ROOT, interpreter


@pytest.mark.integration
@pytest.mark.parametrize(
    "archive",
    sorted((ROOT / "core_genome_archive").glob("*.gb.bz2")),
    ids=lambda p: p.name,
)
def test_reference_genome_regions(archive, tmp_path):
    if os.environ.get("COPRARNA_REFERENCE_TESTS") != "1":
        pytest.skip(
            "Set COPRARNA_REFERENCE_TESTS=1 to run the full bundled-genome oracle"
        )
    perl = interpreter("perl")
    if not perl:
        pytest.skip("Legacy Perl interpreter unavailable")
    genome = tmp_path / archive.name.removesuffix(".bz2")
    genome.write_bytes(bz2.decompress(archive.read_bytes()))
    for region in ("5utr", "3utr", "cds"):
        actual = parse_regions(genome, 200, 100, region).encode()
        result = subprocess.run(
            [
                perl,
                str(ROOT / "coprarna_aux/parse_region_from_genome.pl"),
                str(genome),
                "200",
                "100",
                region,
            ],
            capture_output=True,
            timeout=120,
        )
        assert result.returncode == 0, result.stderr.decode(errors="replace")
        assert actual == result.stdout, f"{archive.name}, {region}"
