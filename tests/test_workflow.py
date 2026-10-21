"""Offline full-workflow regression, including persisted scientific/report data."""

import json
import os
from pathlib import Path
import shutil
import tarfile

import numpy as np
import pytest

from coprarna import alignments, conservation, plots
from coprarna.pipeline import Config, run_prediction
from coprarna.rcompat import load_rdata
from tests.support.legacy import FIXTURES
from tools.build_pipeline_fixture import build


@pytest.fixture
def oracle(tmp_path):
    with tarfile.open(FIXTURES / "workflow/legacy.tar.gz") as archive:
        archive.extractall(tmp_path, filter="data")
    return tmp_path


def test_conservation_and_all_203_jalview_exports_match(oracle, monkeypatch):
    paths = [p for p in (oracle / "evo_alignments").rglob("*") if p.is_file()]
    expected = {p.relative_to(oracle): p.read_bytes() for p in paths}
    matrices = {
        name: load_rdata(oracle / (name + ".Rdata"), name)
        for name in ("conservation_table", "interaction_positions")
    }
    import rdata

    expected_consensus = rdata.read_rda(oracle / "consensus_positions.Rdata")[
        "consensus_both"
    ]
    # Alignment engines are independently exercised in the optional workflow test.
    # This test isolates the Python consensus/mapping algorithm with oracle alignments.
    monkeypatch.setattr(conservation, "command", lambda *args, **kwargs: "")
    conservation.compute_conservation("NC_000913", oracle)
    for path, content in expected.items():
        assert (oracle / path).read_bytes() == content, path
    for name, values in matrices.items():
        for expected_matrix, actual_matrix in zip(
            values, load_rdata(oracle / (name + ".Rdata"), name)
        ):
            np.testing.assert_array_equal(actual_matrix, expected_matrix)

    actual_consensus = rdata.read_rda(oracle / "consensus_positions.Rdata")[
        "consensus_both"
    ]

    def compare(a, b):
        if isinstance(a, dict):
            assert list(a) == list(b)
            for key in a:
                compare(a[key], b[key])
        elif isinstance(a, list):
            assert len(a) == len(b)
            for x, y in zip(a, b):
                compare(x, y)
        else:
            np.testing.assert_array_equal(a, b)

    compare(actual_consensus, expected_consensus)


def test_parallel_alignment_inputs_and_outputs_match(oracle, monkeypatch):
    expected = {
        p.name: p.read_bytes() for p in (oracle / "target_alignments").iterdir()
    }
    shutil.rmtree(oracle / "target_alignments")

    def fake_mafft(args, directory, output):
        source = Path(args[-1])
        assert (directory / source).read_bytes() == expected[source.name]
        (directory / output).write_bytes(expected[Path(output).name])

    monkeypatch.setattr(alignments, "command", fake_mafft)
    alignments.align_targets(
        oracle / "CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv", oracle
    )


def scientific_tools():
    search = (
        os.environ.get("COPRARNA_TOOL_PATH", "")
        + os.pathsep
        + os.environ.get("PATH", "")
    )
    return all(shutil.which(name, path=search) for name in ("mafft", "clustalo"))


@pytest.mark.skipif(not scientific_tools(), reason="MAFFT and Clustal Omega required")
def test_plot_geometry_and_heatmap_match_original_R(oracle):
    data, *_ = plots.region_plot_data(oracle / "CopraRNA_result_all.csv", 20, oracle)
    for feature, points in data.items():
        actual = np.array(
            [[v["start"], v["end"], v["y"], *v["interaction"]] for v in points]
        )
        np.testing.assert_array_equal(
            actual, np.loadtxt(oracle / f"{feature}_geometry.tsv")
        )
        expected = (oracle / f"{feature}_gaps.tsv").read_text().splitlines()
        actual = [",".join(str(int(v)) for v in point["gaps"]) for point in points]
        assert actual == expected
    scores, notes, tree, *_ = plots.evolutionary_heatmap_data(
        oracle, inputfile="CopraRNA2_final_all_ooi.csv"
    )
    np.testing.assert_array_equal(scores, np.loadtxt(oracle / "heatmap_scores.tsv"))
    np.testing.assert_array_equal(
        notes,
        np.loadtxt(oracle / "heatmap_notes.tsv", dtype=str, delimiter="\t", ndmin=2),
    )
    from scipy.cluster.hierarchy import leaves_list

    np.testing.assert_array_equal(
        leaves_list(tree) + 1, np.loadtxt(oracle / "heatmap_order.txt")
    )


@pytest.mark.skipif(
    os.environ.get("COPRARNA_PIPELINE_TESTS") != "1",
    reason="opt-in scientific integration",
)
def test_full_native_prediction_matches_legacy(oracle, tmp_path):
    job = tmp_path / "native"
    build(job)
    expected = json.loads((FIXTURES / "workflow/provenance.json").read_text())
    run_prediction(
        Path("input_sRNA.fa"),
        job,
        Config(winsize=80, maxbpdist=60, cores=2, topcount=20, noclean=True),
    )
    for name in expected["exact_outputs"]:
        assert (job / name).read_bytes() == (oracle / name).read_bytes(), name
    for filename in expected["semantic_Rdata"]:
        name = Path(filename).stem
        for actual, old in zip(
            load_rdata(job / filename, name), load_rdata(oracle / filename, name)
        ):
            np.testing.assert_array_equal(actual, old)


@pytest.mark.skipif(
    os.environ.get("COPRARNA_PIPELINE_TESTS") != "1",
    reason="opt-in scientific integration",
)
def test_full_seeded_copra1_matches_legacy(tmp_path):
    job = tmp_path / "native"
    build(job)
    run_prediction(
        Path("input_sRNA.fa"),
        job,
        Config(
            winsize=80,
            maxbpdist=60,
            cores=2,
            topcount=20,
            cop1=True,
            seed=731,
            noclean=True,
        ),
    )
    with tarfile.open(FIXTURES / "copra1-workflow/legacy.tar.gz") as archive:
        for member in archive.getmembers():
            assert (job / member.name).read_bytes() == archive.extractfile(
                member
            ).read(), member.name
