"""Offline full-workflow regression, including persisted scientific/report data."""

import os
from pathlib import Path
import shutil
import tarfile

import numpy as np
import pytest

from coprarna import alignments, conservation, plots
from coprarna.pipeline import Config, run_prediction
from coprarna.rcompat import load_rdata
from tests.support.fixtures import FIXTURES
from tools.build_pipeline_fixture import build


@pytest.fixture
def oracle(tmp_path):
    with tarfile.open(FIXTURES / "workflow/captured.tar.gz") as archive:
        archive.extractall(tmp_path, filter="data")
    shutil.copyfile(
        FIXTURES / "workflow/organisms.json", tmp_path / "genome_metadata.json"
    )
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
def test_file_and_inline_parameters_produce_identical_predictions(
    tmp_path, monkeypatch
):
    from coprarna import intarna, webserver
    from coprarna.runtime import command

    database = tmp_path / "database"
    database.mkdir()
    jobs = [tmp_path / "file-job", tmp_path / "inline-job"]
    for job in jobs:
        build(job, database)
    for path in jobs[0].glob("*.gb.gz"):
        path.rename(database / path.name)
    for path in jobs[1].glob("*.gb.gz"):
        path.unlink()
    parameter_text = "tAccW=80\ntAccL=60\n"
    supplied = jobs[0] / "parameters.cfg"
    supplied.write_text(parameter_text)
    original = {job: (job / "input_sRNA.fa").read_bytes() for job in jobs}
    calls = []

    def record(arguments, directory, **kwargs):
        if arguments[0] == "IntaRNAsTar":
            path = Path(arguments[arguments.index("--parameterFile") + 1])
            calls.append((arguments[:], path, path.read_text()))
        return command(arguments, directory, **kwargs)

    monkeypatch.setattr(intarna, "command", record)
    monkeypatch.setattr(webserver, "command", record)
    for index, job in enumerate(jobs):
        setting = (
            {"intarna_params": str(supplied)}
            if index == 0
            else {"intarna_params_text": parameter_text}
        )
        run_prediction(
            Path("input_sRNA.fa"),
            job,
            Config(
                cores=2,
                topcount=20,
                noclean=True,
                websrv=True,
                genome_db=str(database),
                **setting,
            ),
        )
        assert (job / "input_sRNA.fa").read_bytes() == original[job]
        assert len(list(job.glob("*.fa.intarna.csv.gz"))) == 3
        assert not list(job.glob("*.fa.intarna.csv"))
        assert not list(job.glob(".intarna-*"))
        assert all(
            Path(line).parent == database
            for line in (job / "genome_paths.txt").read_text().splitlines()
        )
        assert not list(job.glob("*.gb*"))
        for pdf in job.glob("*.pdf"):
            assert pdf.with_suffix(".svg").is_file()
    assert calls and all(contents == parameter_text for _, _, contents in calls)
    assert len(calls) == 8  # Three organisms and one web export in each job.
    assert all(
        path.exists() if path == supplied else not path.exists() for _, path, _ in calls
    )
    for name in (
        "CopraRNA_result.csv",
        "CopraRNA_result_all.csv",
        "coprarna_websrv_table.csv",
    ):
        assert (jobs[0] / name).read_bytes() == (jobs[1] / name).read_bytes(), name


@pytest.mark.skipif(
    os.environ.get("COPRARNA_PIPELINE_TESTS") != "1",
    reason="opt-in scientific integration",
)
def test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors(
    tmp_path,
):
    import csv
    import json
    import zipfile
    import xml.etree.ElementTree as ET

    database = tmp_path / "database"
    database.mkdir()
    job = tmp_path / "native"
    build(job, database, accessions=["NZ_CP999901", "NC_999902", "NC_999903"])
    # Exercise an assembly FASTA ID and a long, previously unknown NZ_ report ID
    # through real native tools, without an availability catalog or network.
    supplied = job / "input_sRNA.fa"
    supplied.write_text(
        supplied.read_text().replace(">NZ_CP999901", ">GCF_000000001.1")
    )
    metadata = json.loads((database / "metadata/NZ_CP999901.json").read_text())
    metadata["query"] = "GCF_000000001.1"
    (database / "metadata/GCF_000000001.1.json").write_text(json.dumps(metadata))
    (job / "cluster.tab").unlink()
    for path in job.glob("*.gb.gz"):
        path.rename(database / path.name)
    original = job / "annotations.fa"
    (job / "input_sRNA.fa").rename(original)
    original_bytes = original.read_bytes()
    run_prediction(original, job, Config(cores=2, websrv=True, genome_db=str(database)))
    assert original.read_bytes() == original_bytes
    with (job / "CopraRNA_result_all.csv").open() as stream:
        rows = list(csv.reader(stream))
    assert len(rows) > 3
    assert rows[0][2:5] == ["NZ_CP999901", "NC_999902", "NC_999903"]
    assert "version:CopraRNA 3.0.0\n" in (job / "CopraRNA_option_file.txt").read_text()
    assert (job / "cluster.tab").read_text().startswith("#id\t")
    with (job / "coprarna_websrv_table.csv").open() as stream:
        web = list(csv.reader(stream))
    assert len(web) == 1 + sum(bool(row[2]) for row in rows[1:])
    assert not (job / "cluster.members").exists()
    assert not (job / "target_alignments").exists()
    assert len(list((job / "IntaRNA").glob("*.csv.gz"))) == 3
    tree = job / "Phylogeny/compatible.treefile"
    assert tree.is_file() and tree.with_suffix(".treefile.svg").is_file()
    for pdf in job.rglob("*.pdf"):
        svg = pdf.with_suffix(".svg")
        assert svg.is_file(), pdf
        root = ET.parse(svg).getroot()
        assert not any(element.tag.endswith("image") for element in root.iter()), svg
    assert not list(job.rglob("*.ps")) and not list(job.rglob("*.png"))
    with zipfile.ZipFile(job / "evo_alignments.zip") as archive:
        assert any(name.endswith("_mapping_result.txt") for name in archive.namelist())


@pytest.mark.skipif(
    os.environ.get("COPRARNA_PIPELINE_TESTS") != "1",
    reason="opt-in scientific integration",
)
def test_native_perfect_helix_is_removed_before_statistics(tmp_path):
    from coprarna.intarna import filter_perfect_helices
    from coprarna.io import read_text
    from coprarna.runtime import command
    from coprarna.tables import INTARNA_HEADER

    (tmp_path / "target.fa").write_text(">self\n" + "C" * 32 + "\n")
    (tmp_path / "query.fa").write_text(">query\n" + "G" * 32 + "\n")
    command(
        [
            "IntaRNAsTar",
            "--target=target.fa",
            "--query=query.fa",
            "--outMode=C",
            "--outCsvCols=" + INTARNA_HEADER.replace(";", ","),
            "--out=interactions.csv.gz",
        ],
        tmp_path,
    )
    output = tmp_path / "interactions.csv.gz"
    assert len(read_text(output).splitlines()) > 1
    report = filter_perfect_helices(output, {"tAccW": "150", "qAccW": "150"})
    assert report["removed_perfect_helices"] >= 1
    assert report["retained_interactions"] == 0
    assert len(read_text(output).splitlines()) == 1
