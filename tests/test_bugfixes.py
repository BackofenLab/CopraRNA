"""Dedicated regressions for each isolated correction on the phase-6 branch."""

from tests.support.fixtures import FIXTURES


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


def test_reverse_flanks_complement_iupac_ambiguity_and_case():
    from coprarna.genome import _reverse_complement

    assert _reverse_complement("ARYKMBVDHN") == "NDHBVKMRYT"
    assert _reverse_complement("arykmbvdhn") == "ndhbvkmryt"


def test_protein_mapping_uses_basename_for_absolute_input(tmp_path):
    from coprarna.genome import extract_proteins

    path = tmp_path / "NC_000001.gb"
    _genbank(path, ["A" * 1000])
    assert extract_proteins(path, organism="abcd") == ">abcd:tag0\nMAA\n"


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


def test_domclust_options_precede_positional_files(tmp_path, monkeypatch):
    from coprarna import homology

    calls = []
    monkeypatch.setattr(
        homology, "command", lambda args, directory, output: calls.append(args)
    )
    (tmp_path / "all.fas.gene").write_text("aaaa tag 80\n")
    (tmp_path / "cluster.members").write_text("1|1 aaaa:tag 1 1 80\n")
    homology.run_domclust(tmp_path)
    assert calls[-1] == [
        "domclust",
        "-HO",
        "-S",
        "-c60",
        "-p0.5",
        "-V0.6",
        "-C80",
        "-o9",
        "all.fas.hom",
        "all.fas.gene",
    ]


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

    return plots, table


def test_lowercase_cds_region_plot_renders(tmp_path, monkeypatch):
    plots, table = _plot_fixture(tmp_path, monkeypatch, "cds", 15)
    assert plots.plot_interaction_regions(table, 15, tmp_path) == ""
    assert (
        (tmp_path / "mRNA_regions_with_histogram.pdf").read_bytes().startswith(b"%PDF")
    )
    assert (tmp_path / "sRNA_regions_with_histogram.svg").is_file()


def test_short_predictions_plot_without_padding_or_out_of_bounds(tmp_path, monkeypatch):

    plots, table = _plot_fixture(tmp_path, monkeypatch, "5utr", 3)
    assert plots.plot_interaction_regions(table, 200, tmp_path) == ""
    assert (
        (tmp_path / "mRNA_regions_extended_list.pdf").read_bytes().startswith(b"%PDF")
    )


def test_domclust_supported_output_preserves_clusters_domains_and_missing_species():
    from coprarna.homology import domclust_table

    genes = "aaaa first 80\nbbbb second 90\ncccc third 70\n"
    members = "1|1 aaaa:first 1 1 40\n1|1 bbbb:second 1 1 90\n2|2 aaaa:first 2 41 80\n2|2 cccc:third 1 1 70\n"
    assert (
        domclust_table(members, genes)
        == "#id\taaaa\tbbbb\tcccc\n1\taaaa:first(1)\tbbbb:second\t\n2\taaaa:first(2)\t\tcccc:third\n"
    )
