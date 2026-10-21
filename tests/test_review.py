"""Explicit contracts for Martin Raden's CopraRNA2 review requests."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import fields
import csv
import gzip
import json
from io import BytesIO
import xml.etree.ElementTree as ET

import pytest

from coprarna import genome, intarna, ncbi, parameters, pipeline, phylogeny
from coprarna.cli import main
from coprarna.helpers import HELPERS
from coprarna.io import open_text, read_text, write_text
from tests.support.legacy import FIXTURES


def test_retired_mode_and_native_parameters_are_not_advertised(capsys):
    names = {field.name for field in fields(pipeline.Config)}
    assert names.isdisjoint({"cop1", "seed", "rcsize", "winsize", "maxbpdist"})
    assert "join-copra1" not in HELPERS and "filter-clusters" not in HELPERS
    with pytest.raises(SystemExit) as result:
        main(["run", "--help"])
    assert result.value.code == 0
    text = capsys.readouterr().out
    for option in ("--cop1", "--seed", "--winsize", "--maxbpdist"):
        assert option not in text
    assert "--intarna-params" in text and "--intarna-params-text" in text


def test_file_parameters_are_forwarded_unchanged(tmp_path):
    path = tmp_path / "parameters.cfg"
    content = "tAccW=80\ntAccL=60\nseedBP=6\n"
    path.write_text(content)
    with parameters.parameter_file(tmp_path, path) as actual:
        assert actual == path.resolve()
        (tmp_path / "CopraRNA_option_file.txt").write_text(
            f"intarna parameter file:{actual}\n"
        )
        assert parameters.arguments(tmp_path) == ["--parameterFile", str(actual)]
        assert parameters.values(tmp_path)["tAccW"] == "80"
    assert path.read_text() == content
    assert (tmp_path / "IntaRNA_parameters.txt").read_text() == content


@pytest.mark.parametrize("fail", [False, True])
def test_inline_parameters_are_temporary_even_after_failure(tmp_path, fail):
    content = "tAccW=80\nqAccW=0\n"
    with (
        pytest.raises(RuntimeError) if fail else __import__("contextlib").nullcontext()
    ):
        with parameters.parameter_file(tmp_path, text=content) as path:
            assert path.read_text() == content
            if fail:
                raise RuntimeError("prediction failed")
    assert not path.exists()
    assert (tmp_path / "IntaRNA_parameters.txt").read_text() == content


def test_parameter_inputs_are_mutually_exclusive(tmp_path):
    with pytest.raises(ValueError, match="not both"):
        with parameters.parameter_file(tmp_path, "a.cfg", "tAccW=80"):
            pass
    with pytest.raises(SystemExit):
        main(["run", "--intarna-params", "a.cfg", "--intarna-params-text", "tAccW=80"])


def test_intarna_older_than_three_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(parameters, "command", lambda *args: "IntaRNA 2.2.0\n")
    with pytest.raises(RuntimeError, match="3 or newer"):
        parameters.check_version(tmp_path)


def row(pairs, query=100, target=300, structure=None):
    return {
        "seq1": "C" * target,
        "seq2": "G" * query,
        "hybridDP": structure or "(" * pairs + "&" + ")" * pairs,
    }


@pytest.mark.parametrize(
    "record,windows,expected",
    [
        (row(90), {}, True),
        (row(89), {}, False),
        (row(100), {}, True),
        (row(90, query=300, target=100), {}, True),
        (row(72, query=300), {"tAccW": "80"}, True),
        (row(72, query=300), {"tAccW": "0"}, False),
        (row(100, structure="(" * 50 + "." + "(" * 50 + "&" + ")" * 100), {}, False),
        (row(100, structure="(" * 100 + "&" + ")" * 40 + "." + ")" * 60), {}, False),
        (row(100, structure="(" * 100 + "&" + ")" * 99), {}, False),
    ],
)
def test_full_length_stacking_filter(record, windows, expected):
    assert intarna.perfect_helix(record, windows) is expected


def test_gzip_artifacts_are_removed_before_gev_fit(tmp_path, monkeypatch):
    from coprarna import energy, tables

    source = tmp_path / "sample.fa.intarna.csv.gz"
    header = tables.INTARNA_HEADER.split(";")
    records = []
    for tag, pairs, score in [
        ("self", 100, -1000),
        ("other", 50, -12),
        ("other2", 40, -10),
        ("other3", 45, -11),
    ]:
        record = dict.fromkeys(header, "0")
        record.update(
            row(pairs),
            id1=tag,
            id2="query",
            E=str(score),
            E_norm=str(score),
            start1="1",
            end1=str(pairs),
            start2="1",
            end2=str(pairs),
        )
        records.append(record)
    with open_text(source, "wt") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=header, delimiter=";", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)
    assert intarna.filter_perfect_helices(source) == {
        "removed_perfect_helices": 1,
        "retained_interactions": 3,
        "removed_interactions": [
            {
                "id1": "self",
                "id2": "query",
                "start1": "1",
                "end1": "100",
                "start2": "1",
                "end2": "100",
            }
        ],
    }
    write_text(
        tmp_path / "sample.fa.intarna.sorted.csv.gz",
        "".join(intarna.numeric_sort(read_text(source).splitlines(keepends=True), 14)),
    )
    tables.disentangle(tmp_path)
    (tmp_path / "CopraRNA_option_file.txt").write_text("region:5utr\n")
    fitted = []
    monkeypatch.setattr(
        energy, "fit_gev", lambda values: fitted.append(values) or [0.1, 1, 0]
    )
    energy.add_energy_pvalues(tmp_path)
    assert fitted == [[12.0, 11.0, 10.0]]
    assert "self" not in read_text(source)
    assert not list(tmp_path.glob(".filtered-*"))


def test_compressed_genomes_parse_without_a_plain_copy(tmp_path):
    source = FIXTURES / "region-linear-5utr/input/NC_000001.gb"
    compressed = tmp_path / "NC_000001.gb.gz"
    compressed.write_bytes(gzip.compress(source.read_bytes()))
    for region in ("5utr", "3utr", "cds"):
        assert genome.parse_regions(
            compressed, 200, 100, region
        ) == genome.parse_regions(source, 200, 100, region)
    assert not (tmp_path / "NC_000001.gb").exists()
    assert genome.check_gene_cds(compressed) == genome.check_gene_cds(source)


def test_shared_genome_cache_is_read_directly(tmp_path, monkeypatch):
    database = tmp_path / "db"
    database.mkdir()
    job = tmp_path / "job"
    job.mkdir()
    cached = database / "NC_000001.gb.gz"
    cached.write_bytes(gzip.compress(b"genome\n//\n"))
    monkeypatch.setattr(
        ncbi, "download_genbank", lambda *args: pytest.fail("cached genome downloaded")
    )
    assert ncbi.genome_path("NC_000001", job, database) == cached
    assert not list(job.iterdir())


def test_parallel_jobs_download_a_shared_genome_once(tmp_path, monkeypatch):
    database = tmp_path / "db"
    calls = []

    def download(identifier, path):
        calls.append(identifier)
        path.write_bytes(gzip.compress(b"genome\n//\n"))

    monkeypatch.setattr(ncbi, "download_genbank", download)
    with ThreadPoolExecutor(max_workers=4) as pool:
        paths = list(
            pool.map(
                lambda _: ncbi.genome_path("NC_000001", tmp_path, database), range(4)
            )
        )
    assert calls == ["NC_000001"]
    assert len(set(paths)) == 1


@pytest.mark.parametrize("server_gzip", [False, True])
def test_download_is_streamed_and_kept_compressed(tmp_path, server_gzip):
    content = b"LOCUS test\n" + b"ACGT" * 500000 + b"\n//\n"

    class Stream(BytesIO):
        def read(self, size=-1):
            assert size >= 0, "download buffered the entire genome"
            return super().read(size)

    stream = Stream(gzip.compress(content) if server_gzip else content)
    stream.headers = {"Content-Encoding": "gzip"} if server_gzip else {}
    responses = [
        BytesIO(b"<eSearchResult><IdList><Id>1</Id></IdList></eSearchResult>"),
        BytesIO(b"<eSummaryResult/>"),
        stream,
    ]

    def opener(request):
        if hasattr(request, "headers"):
            assert request.get_header("Accept-encoding") == "gzip"
        return responses.pop(0)

    target = tmp_path / "NC_000001.gb.gz"
    ncbi.download_genbank("NC_000001", target, opener=opener)
    assert gzip.decompress(target.read_bytes()) == content
    assert not list(tmp_path.glob("*.gb"))


def test_truncated_download_is_not_cached(tmp_path):
    responses = [
        BytesIO(b"<eSearchResult><IdList><Id>1</Id></IdList></eSearchResult>"),
        BytesIO(b"<eSummaryResult/>"),
        BytesIO(b"LOCUS unfinished\n"),
    ]
    target = tmp_path / "NC_000001.gb.gz"
    with pytest.raises(ValueError, match="Incomplete"):
        ncbi.download_genbank("NC_000001", target, opener=lambda _: responses.pop(0))
    assert not target.exists() and not list(tmp_path.glob(".NC_*"))


@pytest.mark.parametrize("name", ["input_sRNA.fa", "annotations.fa"])
def test_original_fasta_is_not_modified_or_moved(tmp_path, monkeypatch, name):
    content = b">NC_000913 description-with-hyphen\nacguacgu\n>NC_000964\nacguacgu\n>NC_002516\nacguacgu\n"
    original = tmp_path / name
    original.write_bytes(content)
    monkeypatch.setattr(parameters, "check_version", lambda _: "IntaRNA 3.4.1")

    def stop(*args, **kwargs):
        raise RuntimeError("stop after preparation")

    monkeypatch.setattr(pipeline.homology, "run_homology", stop)
    with pytest.raises(RuntimeError, match="after preparation"):
        pipeline.run_prediction(
            original, tmp_path, pipeline.Config(intarna_params_text="tAccW=80\n")
        )
    assert original.read_bytes() == content
    assert not list(tmp_path.glob(".intarna-*"))
    assert list((tmp_path / "FASTA").glob("coprarna_input_*.fa"))


def test_tree_svg_uses_only_python_and_vector_paths(tmp_path):
    tree = tmp_path / "compatible.treefile"
    tree.write_text("(NC_000001:0.1,(NC_000002:0.2,NC_000003:0.3):0.4);\n")
    output = tmp_path / "compatible.treefile.svg"
    phylogeny.render_tree(tree, output)
    root = ET.parse(output).getroot()
    assert root.tag.endswith("svg")
    assert any(node.tag.endswith("path") for node in root.iter())
    assert not any(node.tag.endswith("image") for node in root.iter())


@pytest.mark.parametrize("name", ["parameters.cfg", "IntaRNA_parameters.txt"])
def test_parameter_snapshot_preserves_original_bytes(tmp_path, name):
    original = tmp_path / name
    content = b"tAccW=80\r\ntAccL=60\r\n"
    original.write_bytes(content)
    with parameters.parameter_file(tmp_path, original):
        assert (tmp_path / "IntaRNA_parameters.txt").read_bytes() == content
    assert original.read_bytes() == content


def test_prediction_cache_signature_tracks_contents_not_paths(tmp_path):
    (tmp_path / "target.fa").write_text(">target\nACGU\n")
    (tmp_path / "query.fa").write_text(">query\nUGCA\n")
    (tmp_path / "CopraRNA_option_file.txt").write_text(
        "IntaRNA version:IntaRNA 3.4.1\n"
    )
    first = intarna.prediction_signature(tmp_path, "target.fa", "query.fa")
    other = tmp_path / "other"
    import shutil

    shutil.copytree(tmp_path, other)
    assert intarna.prediction_signature(other, "target.fa", "query.fa") == first
    (tmp_path / "query.fa").write_text(">query\nUGCC\n")
    assert intarna.prediction_signature(tmp_path, "target.fa", "query.fa") != first
    (tmp_path / "query.fa").write_text(">query\nUGCA\n")
    (tmp_path / "params.cfg").write_text("tAccW=80\n")
    (tmp_path / "CopraRNA_option_file.txt").write_text(
        f"IntaRNA version:IntaRNA 3.4.1\nintarna parameter file:{tmp_path / 'params.cfg'}\n"
    )
    assert intarna.prediction_signature(tmp_path, "target.fa", "query.fa") != first


def test_report_retains_missing_comparative_pvalues(tmp_path):
    from coprarna.webserver import prepare_webserver_output
    import shutil

    fixture = FIXTURES / "webserver-output"
    shutil.copytree(fixture / "input", tmp_path, dirs_exist_ok=True)
    spec = json.loads((fixture / "case.json").read_text())
    result, interactions = [tmp_path / value for value in spec["args"][:2]]
    lines = result.read_text().splitlines()
    fields = lines[1].split(",")
    fields[:2] = ["NA", "NA"]
    lines[1] = ",".join(fields)
    result.write_text("\n".join(lines) + "\n")
    prepare_webserver_output(result, interactions, tmp_path)
    assert (tmp_path / "coprarna_internal_table.csv").read_text().splitlines()[1].split(
        ","
    )[1:3] == ["NA", "NA"]


def test_all_extended_plot_pages_are_exported_as_vectors(tmp_path, monkeypatch):
    from tests.test_bugfixes import _plot_fixture

    plots, table = _plot_fixture(tmp_path, monkeypatch, "5utr", 30)
    plots.plot_interaction_regions(table, 30, tmp_path)
    for name in ("mRNA_regions_extended_list", "sRNA_regions_with_extended_list"):
        drawing = ET.parse(tmp_path / (name + ".svg")).getroot()
        assert sum(node.tag.endswith("svg") for node in drawing.iter()) == 3
        assert not any(node.tag.endswith("image") for node in drawing.iter())
        pdf = (tmp_path / (name + ".pdf")).read_bytes()
        assert b"/Count 2" in pdf


def test_functional_enrichment_has_vector_pdf_and_svg(tmp_path):
    from coprarna.plots import plot_functional_enrichment

    prefix = tmp_path / "enriched_heatmap_big"
    plot_functional_enrichment(FIXTURES / "heatmap-json/expected/copraRNA.json", prefix)
    assert prefix.with_suffix(".pdf").read_bytes().startswith(b"%PDF")
    drawing = ET.parse(prefix.with_suffix(".svg")).getroot()
    assert not any(node.tag.endswith("image") for node in drawing.iter())
