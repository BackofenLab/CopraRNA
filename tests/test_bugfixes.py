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
