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


def _genbank(path, sequences):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    records = []
    for i, sequence in enumerate(sequences):
        record = SeqRecord(
            Seq(sequence), id=f"NC_00000{i + 1}.1", name=f"NC_00000{i + 1}"
        )
        record.annotations = {"molecule_type": "DNA", "topology": "linear"}
        record.features = [
            SeqFeature(
                SimpleLocation(400, 550, strand=1),
                type="CDS",
                qualifiers={"locus_tag": [f"tag{i}"], "translation": ["MAA"]},
            )
        ]
        records.append(record)
    SeqIO.write(records, path, "genbank")


def test_region_features_use_their_own_genbank_record(tmp_path):
    from coprarna.genome import parse_regions

    path = tmp_path / "records.gb"
    _genbank(path, ["A" * 1000, "C" * 1000])
    actual = parse_regions(path, 100, 50, "5utr")
    assert actual == ">tag0\n" + "A" * 150 + "\n>tag1\n" + "C" * 150 + "\n"
    with legacy_behavior():
        assert parse_regions(path, 100, 50, "5utr").endswith("A" * 150 + "\n")


def test_reverse_flanks_complement_iupac_ambiguity_and_case():
    from coprarna.genome import _reverse_complement

    assert _reverse_complement("ARYKMBVDHN") == "NDHBVKMRYT"
    assert _reverse_complement("arykmbvdhn") == "ndhbvkmryt"
    with legacy_behavior():
        assert _reverse_complement("ARYKMBVDHN") == "NHDVBMKYRT"


def test_protein_mapping_uses_basename_for_absolute_input(tmp_path):
    from coprarna.genome import extract_proteins

    path = tmp_path / "NC_000001.gb"
    _genbank(path, ["A" * 1000])
    assert extract_proteins(path, "abcd\tNC_000001\n") == ">abcd:tag0\nMAA\n"
    with legacy_behavior():
        assert extract_proteins(path, "abcd\tNC_000001\n") == ">:tag0\nMAA\n"


def test_suboptimal_keeps_second_best_energy_in_sorted_input(tmp_path):
    from coprarna.tables import disentangle

    (tmp_path / "sample.fa.intarna.sorted.csv").write_text(
        "id;E\nx;-10\nx;-9\nx;-8\ny;-5\ny;-4\n"
    )
    disentangle(tmp_path)
    assert (tmp_path / "sample_subopt.intarna.csv").read_text().splitlines()[1:] == [
        "x;-9",
        "y;-4",
    ]
    with legacy_behavior():
        disentangle(tmp_path)
    assert (tmp_path / "sample_subopt.intarna.csv").read_text().splitlines()[1:] == [
        "x;-8",
        "y;-4",
    ]


def test_heatmap_json_escapes_quotes_and_backslashes(tmp_path):
    import json
    from coprarna.exports import heatmap_json

    source = next((FIXTURES / "heatmap-json/input").iterdir())
    rows = [line.split("\t") for line in source.read_text().splitlines()]
    rows[2][11] = 'term "quoted" \\ slash'
    rows[3][4] = 'gene "quoted" \\ slash'
    path = tmp_path / "enrichment.txt"
    path.write_text("\n".join("\t".join(row) for row in rows) + "\n")
    heatmap_json(path)
    data = json.loads((tmp_path / "copraRNA.json").read_text())
    assert data[0]["feature"] == rows[2][11]
    assert data[0]["source"] == rows[3][4]


def test_copra1_rho_error_uses_histogram_density():
    import numpy as np
    import pytest
    from coprarna.combination import _histogram_error

    # Independently evaluated by R hist(...,breaks=100): sum((1-h$density)^2).
    skew = np.r_[np.zeros(90), np.linspace(0.1, 1, 10)]
    assert _histogram_error(skew) == pytest.approx(8010.0, abs=1e-9)
    assert _histogram_error(np.linspace(0.005, 0.995, 100)) < 1e-20
    with legacy_behavior():
        assert _histogram_error(skew) == 0


def test_domclust_options_precede_positional_files(tmp_path, monkeypatch):
    from coprarna import homology

    calls = []
    monkeypatch.setattr(
        homology, "command", lambda args, directory, output: calls.append(args)
    )
    homology.run_domclust(tmp_path)
    assert calls[-1] == [
        "domclust",
        "-HO",
        "-S",
        "-c60",
        "-p0.5",
        "-V0.6",
        "-C80",
        "-o5",
        "all.fas.hom",
        "all.fas.gene",
    ]
    with legacy_behavior():
        homology.run_domclust(tmp_path)
    assert calls[-1][1:3] == ["all.fas.hom", "all.fas.gene"]


def _plot_fixture(tmp_path, monkeypatch, region, count):
    from coprarna import plots

    (tmp_path / "CopraRNA_option_file.txt").write_text(
        f"nt upstream:200\nnt downstream:100\nregion:{region}\nwebsrv:0\n"
    )
    table = tmp_path / "predictions.csv"
    table.write_text(
        "fdr,p,org,Annotation\n"
        + "".join("0,0,tag(g|-1|.01|2|4|3|5),annotation\n" for _ in range(count))
    )

    def geometry(path, number, directory):
        points = [
            dict(
                start=40.0,
                end=340.0,
                interaction=[200.0, 210.0],
                gaps=[],
                y=1 + i * 2,
                cluster=i,
            )
            for i in range(number)
        ]
        return (
            {"mRNA": points, "sRNA": points},
            [str(i) for i in range(number)],
            ["org"],
        )

    monkeypatch.setattr(plots, "region_plot_data", geometry)

    def ps_output(args, directory):
        output = next(
            value.split("=", 1)[1]
            for value in args
            if value.startswith("-sOutputFile=")
        )
        (directory / output).write_text("%!PS-Adobe-3.0\n")

    monkeypatch.setattr(plots, "command", ps_output)
    return plots, table


def test_lowercase_cds_region_plot_renders(tmp_path, monkeypatch):
    plots, table = _plot_fixture(tmp_path, monkeypatch, "cds", 15)
    assert plots.plot_interaction_regions(table, 15, tmp_path) == ""
    assert (
        (tmp_path / "mRNA_regions_with_histogram.pdf").read_bytes().startswith(b"%PDF")
    )
    assert (
        (tmp_path / "sRNA_regions_with_histogram.png")
        .read_bytes()
        .startswith(b"\x89PNG")
    )


def test_short_predictions_plot_without_padding_or_out_of_bounds(tmp_path, monkeypatch):
    import pytest

    plots, table = _plot_fixture(tmp_path, monkeypatch, "5utr", 3)
    assert plots.plot_interaction_regions(table, 200, tmp_path) == ""
    assert (
        (tmp_path / "mRNA_regions_extended_list.pdf").read_bytes().startswith(b"%PDF")
    )
    with legacy_behavior(), pytest.raises(ValueError, match="beyond available"):
        plots.plot_interaction_regions(table, 200, tmp_path)
