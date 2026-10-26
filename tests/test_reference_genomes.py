"""Identity, provenance and extraction contracts for refreshed NCBI references."""

import hashlib
from io import StringIO
import json
import os

from Bio import SeqIO
import pytest

from coprarna import genome, organisms
from coprarna.io import genbank_records
from tests.support.fixtures import ROOT

ARCHIVE = ROOT / "core_genome_archive"


def entries():
    return json.loads((ARCHIVE / "manifest.json").read_text())["genomes"]


def test_reference_archive_has_17_current_compressed_genomes():
    records = entries()
    assert len(records) == 17
    assert {p.name for p in ARCHIVE.glob("*.gb.gz")} == {r["file"] for r in records}
    assert not list(ARCHIVE.glob("*.bz2"))
    for item in records:
        assert (
            hashlib.sha256((ARCHIVE / item["file"]).read_bytes()).hexdigest()
            == item["sha256"]
        )
        metadata = organisms.resolve(item["query"], ARCHIVE)
        assert item["accession"] in metadata["replicons"]
        assert item["assembly"] == metadata["assembly"]


@pytest.mark.integration
@pytest.mark.skipif(
    os.environ.get("COPRARNA_REFERENCE_TESTS") != "1",
    reason="opt-in reference extraction",
)
@pytest.mark.parametrize("path", sorted(ARCHIVE.glob("*.gb.gz")), ids=lambda p: p.name)
def test_reference_genome_regions(path):
    (record,) = genbank_records(path)
    item = next(r for r in entries() if r["file"] == path.name)
    assert record.id == item["accession"] and len(record) == item["length"]
    assert record.annotations["date"] == item["record_date"]
    cds = [f for f in record.features if f.type == "CDS"]
    assert len(cds) == item["cds_count"] > 0
    expected = [
        (f.qualifiers["locus_tag"][0], f.qualifiers["translation"][0])
        for f in cds
        if "translation" in f.qualifiers and "locus_tag" in f.qualifiers
    ]
    actual = list(SeqIO.parse(StringIO(genome.extract_proteins(path)), "fasta"))
    assert [(r.id.split(":")[1], str(r.seq)) for r in actual] == expected
    tags = {f.qualifiers.get("locus_tag", [""])[0] for f in cds}
    for region in ("5utr", "3utr", "cds"):
        parsed = list(
            SeqIO.parse(StringIO(genome.parse_regions(path, 200, 100, region)), "fasta")
        )
        assert parsed
        assert all(r.id in tags and len(r) > 0 for r in parsed)
        assert all(set(str(r.seq).upper()) <= set("ACGTRYSWKMBDHVN") for r in parsed)
        extended = list(genome.extract_targets(path, 200, 100, region, 100))
        assert [(r.tag, r.requested) for r in extended] == [
            (r.id, str(r.seq)) for r in parsed
        ]
        assert all(
            0 <= r.start <= 100 and 0 <= len(r.sequence) - r.end <= 100
            for r in extended
        )
        if record.annotations.get("topology") == "circular":
            assert all(r.start == len(r.sequence) - r.end == 100 for r in extended)
