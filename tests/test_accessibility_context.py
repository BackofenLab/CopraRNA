"""Requested-region contracts independent of the unchanged scientific oracles."""

import csv
from dataclasses import asdict
from io import StringIO
import json
import math
import os
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import pytest

from coprarna import genome, intarna, parameters, pipeline, tables
from coprarna.io import ensure_input_is_distinct, open_text, read_text
from coprarna.runtime import command


def settings(directory, contents=""):
    path = directory / "native.cfg"
    path.write_text(contents)
    (directory / "CopraRNA_option_file.txt").write_text(
        f"intarna parameter file:{path}\ncore count:1\n"
    )
    return path


@pytest.mark.parametrize(
    "contents, expected",
    [
        ("", 100),
        ("tAccL=60\n", 60),
        ("accL=70\n", 70),
        ("tAccL=0\n", 150),
        ("tAccL=0\ntAccW=80\n", 80),
        ("accL=0\naccW=90\n", 90),
        ("accL=0\naccW=0\n", None),
        ("qAccL=30\n", 100),
        ("tAccL=60 # comment\ntAccW=80\n", 60),
    ],
)
def test_context_uses_native_span_semantics(tmp_path, contents, expected):
    path = settings(tmp_path, contents)
    assert parameters.target_context(tmp_path) == expected
    assert path.read_text() == contents


@pytest.mark.parametrize("contents", ["tAccL=-1\n", "accL=60\ntAccL=50\n"])
def test_invalid_or_conflicting_context_fails(tmp_path, contents):
    settings(tmp_path, contents)
    with pytest.raises(ValueError):
        parameters.target_context(tmp_path)


def genbank(directory, *, strand=1, circular=False, start=300, end=500, split=False):
    rng = random.Random(42)
    sequence = "".join(rng.choices("ACGT", k=1000))
    record = SeqRecord(Seq(sequence), id="NC_000001.1", name="NC_000001")
    record.annotations = {
        "molecule_type": "DNA",
        "topology": "circular" if circular else "linear",
    }
    location = SimpleLocation(start, end, strand=strand)
    if split:
        parts = [
            SimpleLocation(start, start + 30, strand=strand),
            SimpleLocation(end - 40, end, strand=strand),
        ]
        location = CompoundLocation(parts if strand == 1 else parts[::-1])
    record.features = [
        SeqFeature(location, type="CDS", qualifiers={"locus_tag": ["gene"]})
    ]
    path = directory / "NC_000001.gb.gz"
    with open_text(path, "wt") as stream:
        SeqIO.write(record, stream, "genbank")
    return path, sequence, location


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("region", ["5utr", "3utr", "cds"])
@pytest.mark.parametrize(
    "start,end", [(300, 500), (20, 220), (800, 980), (0, 200), (800, 1000)]
)
def test_genomic_context_strand_boundaries_and_requested_region(
    tmp_path, strand, circular, region, start, end
):
    path, source, location = genbank(
        tmp_path, strand=strand, circular=circular, start=start, end=end
    )
    up, down, padding = 100, 150, 60
    # Independent coordinates in the oriented genomic sequence. The UTR modes
    # contain exactly up+down bases, with the same boundary as the public helper.
    oriented = source if strand == 1 else str(Seq(source).reverse_complement())
    gene_start = start if strand == 1 else len(source) - end
    gene_end = end if strand == 1 else len(source) - start
    anchor = gene_end if region == "3utr" else gene_start
    left = anchor - up
    right = (gene_end if region == "cds" else anchor) + down

    def extract(a, b):
        if circular:
            return "".join(oriented[i % len(oriented)] for i in range(a, b))
        return oriented[max(0, a) : min(len(oriented), b)]

    requested = extract(left, right)
    records = list(genome.extract_targets(path, up, down, region, padding))
    if not circular and region != "cds" and len(requested) < 150:
        assert records == []
        return
    (actual,) = records
    assert actual.requested == requested
    assert actual.sequence == extract(left - padding, right + padding)
    assert actual.start == (padding if circular else min(padding, max(0, left)))
    fasta = list(
        SeqIO.parse(StringIO(genome.parse_regions(path, up, down, region)), "fasta")
    )
    assert str(fasta[0].seq) == requested


@pytest.mark.parametrize("strand", [1, -1])
def test_spliced_cds_keeps_exons_and_adds_external_context(tmp_path, strand):
    path, source, location = genbank(tmp_path, strand=strand, split=True)
    (actual,) = genome.extract_targets(path, 10, 20, "cds", 40)
    if strand == 1:
        before, after = source[250:300], source[500:560]
    else:
        before = str(Seq(source[500:550]).reverse_complement())
        after = str(Seq(source[240:300]).reverse_complement())
    assert actual.sequence == before + str(location.extract(Seq(source))) + after
    assert (actual.start, actual.end) == (40, 140)


@pytest.mark.parametrize("circular", [False, True])
def test_unlimited_span_and_multiple_circular_wraps(tmp_path, circular):
    path, source, _ = genbank(tmp_path, circular=circular)
    (actual,) = genome.extract_targets(path, 100, 150, "5utr", None)
    if circular:
        assert actual.start == 1000
        assert actual.sequence == "".join(source[i % 1000] for i in range(-800, 1450))
    else:
        assert actual.sequence == source
        assert (actual.start, actual.end) == (200, 450)


def test_context_and_interval_invalidate_cache(tmp_path):
    settings(tmp_path)
    (tmp_path / "targets.fa").write_text(">gene\nCCCCC\n")
    (tmp_path / "query.fa").write_text(">query\nGGGGG\n")
    context = tmp_path / "targets.fa.context.json"
    signatures = []
    for sequence, start, end in [
        ("AAAAACCCCCAAAAA", 5, 10),
        ("UUUUUCCCCCAAAAA", 5, 10),
        ("UUUUUCCCCCAAAAA", 4, 9),
    ]:
        context.write_text(
            json.dumps([asdict(genome.TargetRegion("gene", sequence, start, end))])
        )
        signatures.append(
            intarna.prediction_signature(
                tmp_path,
                "targets.fa",
                "query.fa",
                native_version="IntaRNA 3.4.1",
                context=context,
            )
        )
    assert len(set(signatures)) == 3


def test_context_sidecar_is_protected_and_cleanup_preserves_nested_input(tmp_path):
    name = "NC_000001_upfromstartpos_200_down_100.fa.context.json"
    original = tmp_path / "IntaRNA" / name
    original.parent.mkdir()
    original.write_text("original input\n")
    generated = tmp_path / name
    generated.write_text("generated context\n")
    with pytest.raises(ValueError, match="collides"):
        ensure_input_is_distinct(generated, tmp_path)
    (tmp_path / "CopraRNA_option_file.txt").write_text(f"sRNA FASTA:{original}\n")
    pipeline._cleanup(tmp_path, pipeline.Config())
    assert original.read_text() == "original input\n"
    assert generated.read_text() == "generated context\n"


def test_failed_native_group_preserves_prior_csv_and_removes_temporaries(tmp_path):
    settings(tmp_path)
    output = tmp_path / "result.csv.gz"
    output.write_bytes(b"prior result")

    def fail(*args):
        raise RuntimeError("native failure")

    with pytest.raises(RuntimeError, match="native failure"):
        intarna._predict_targets(
            [genome.TargetRegion("gene", "AAAAACCCCCAAAAA", 5, 10)],
            "query.fa",
            output,
            tmp_path,
            tables.INTARNA_HEADER.replace(";", ","),
            runner=fail,
        )
    assert output.read_bytes() == b"prior result"
    assert not list(tmp_path.glob(".intarna-context-*"))


@pytest.mark.skipif(
    os.environ.get("COPRARNA_PIPELINE_TESTS") != "1", reason="opt-in native integration"
)
def test_native_context_changes_accessibility_and_preserves_report_coordinates(
    tmp_path,
):
    settings(
        tmp_path,
        "tAccW=80\ntAccL=60\nseedBP=3\noutMaxE=100\noutOverlap=B\noutNumber=3\n"
        "tIdxPos0=500\nqIdxPos0=200\ntRegion=500-501\n",
    )
    parameters.check_version(tmp_path)
    # The C-rich requested target can fold with the added G-rich genomic flank.
    target = genome.TargetRegion("gene", "G" * 6 + "AAAA" + "C" * 12 + "A" * 16, 10, 22)
    short = genome.TargetRegion("short", "AAAA" + "C" * 8 + "AAAA", 4, 12)
    (tmp_path / "query.fa").write_text(">query\nGGGGGGGGGGGG\n")
    columns = tables.INTARNA_HEADER.replace(";", ",") + ",bpList"
    calls = []

    def record(args, directory):
        calls.append(args)
        return command(args, directory)

    intarna._predict_targets(
        [target, short], "query.fa", "padded.csv.gz", tmp_path, columns, runner=record
    )
    assert len(calls) == 2
    assert "--tRegion=1-12" in calls[0] and "--tIdxPos0=-10" in calls[0]
    assert "--tRegion=1-8" in calls[1] and "--tIdxPos0=-4" in calls[1]
    rows = list(
        csv.DictReader(StringIO(read_text(tmp_path / "padded.csv.gz")), delimiter=";")
    )
    assert {r["id1"] for r in rows} == {"gene", "short"}
    for row in rows:
        size = 12 if row["id1"] == "gene" else 8
        assert row["seq1"] == "C" * size
        assert 1 <= int(row["start1"]) <= int(row["end1"]) <= size
        for key in ("seedStart1", "seedEnd1"):
            assert all(1 <= int(v) <= size for v in row[key].split(":"))
        assert row["subseqDB"].startswith(row["start1"] + "C")
        assert row["hybridDB"].startswith(row["start1"] + "|")
        assert float(row["E_norm"]) == pytest.approx(
            float(row["E"]) / math.log(size * 12), rel=1e-5
        )
        assert all(
            1 <= int(pair.split(",")[0][1:]) <= size
            for pair in row["bpList"].split(":")
        )
    intarna._predict_targets(
        [genome.TargetRegion("gene", "C" * 12, 0, 12)],
        "query.fa",
        "bare.csv.gz",
        tmp_path,
        columns,
    )
    bare = list(
        csv.DictReader(StringIO(read_text(tmp_path / "bare.csv.gz")), delimiter=";")
    )
    assert all(float(r["ED1"]) == 0 for r in bare)
    assert any(float(r["ED1"]) > 0 for r in rows if r["id1"] == "gene")
    assert not list(tmp_path.glob(".intarna-context-*"))
