"""NCBI resolution contracts without live network access."""

from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
import json
from urllib.parse import parse_qs, urlparse

import pytest

from coprarna import organisms, ncbi
from coprarna.cli import main


def metadata(query="NC_000001", assembly="GCF_000000001.1", replicons=None):
    return dict(
        schema_version=1,
        query=query,
        assembly=assembly,
        organism="Test organism strain",
        replicons=replicons or ["NC_000001.2"],
        retrieved_at="2026-10-22T00:00:00+00:00",
    )


def response(value):
    return BytesIO(json.dumps(value).encode())


def service(paged=False, ambiguous=False):
    urls = []

    def opener(url):
        urls.append(url)
        endpoint = urlparse(url).path
        params = parse_qs(urlparse(url).query)
        if endpoint.endswith("esummary.fcgi"):
            return response(
                {"result": {"uids": ["1"], "1": {"accessionversion": "NC_000001.2"}}}
            )
        if endpoint.endswith("sequence_assemblies"):
            assert "/NC_000001.2/" in endpoint
            return response({"accessions": ["GCF_000000001.1", "GCA_000000001.1"]})
        if endpoint.endswith("dataset_report"):
            reports = [
                dict(
                    accession="GCF_000000001.1",
                    organism={"organism_name": "Test organism strain"},
                    assembly_info={
                        "assembly_level": "Complete Genome",
                        "assembly_status": "current",
                    },
                )
            ]
            if ambiguous:
                reports.append(dict(reports[0], accession="GCF_000000002.1"))
            return response({"reports": reports, "total_count": len(reports)})
        assert endpoint.endswith("sequence_reports")
        sequences = [
            dict(
                refseq_accession="NC_000003.1",
                assigned_molecule_location_type="Plasmid",
                length=100,
                role="assembled-molecule",
            ),
            dict(
                refseq_accession="NC_000001.2",
                assigned_molecule_location_type="Chromosome",
                length=1000,
                role="assembled-molecule",
            ),
            dict(
                refseq_accession="NC_000002.1",
                assigned_molecule_location_type="Chromosome",
                length=500,
                role="assembled-molecule",
            ),
        ]
        if paged:
            if not params.get("page_token"):
                return response(
                    {
                        "reports": sequences[:1],
                        "next_page_token": "two",
                        "total_count": 3,
                    }
                )
            assert params["page_token"] == ["two"]
            return response({"reports": sequences[1:], "total_count": 3})
        return response({"reports": sequences, "total_count": 3})

    return opener, urls


@pytest.mark.parametrize("paged", [False, True])
def test_resolve_complete_assembly_chromosomes_first_then_plasmids_and_cache(
    tmp_path, paged
):
    opener, urls = service(paged=paged)
    actual = organisms.resolve("NC_000001", tmp_path, opener=opener)
    assert actual["replicons"] == ["NC_000001.2", "NC_000002.1", "NC_000003.1"]
    assert actual["organism"] == "Test organism strain"
    assert actual["assembly"] == "GCF_000000001.1"
    assert len(urls) == (5 if paged else 4)
    assert (
        organisms.resolve(
            "NC_000001", tmp_path, opener=lambda _: pytest.fail("cache miss")
        )
        == actual
    )
    assert not list(tmp_path.glob("*.gb*"))


def test_assembly_input_needs_no_sequence_lookup(tmp_path):
    opener, urls = service()
    result = organisms.resolve("GCF_000000001.1", tmp_path, opener=opener)
    assert result["query"] == "GCF_000000001.1"
    assert len(urls) == 2


@pytest.mark.parametrize(
    "identifier",
    [
        "NC_000001extra",
        "../NC_000001",
        "NZ_",
        "GCA_000000001.1",
        "NC_001.0",
        "NC_001\n",
        "bad",
    ],
)
def test_bad_accession_rejected_before_side_effects(tmp_path, identifier):
    with pytest.raises(ValueError, match="Invalid RefSeq"):
        organisms.resolve(
            identifier, tmp_path / "db", opener=lambda _: pytest.fail("network request")
        )
    assert not (tmp_path / "db").exists()


def test_ambiguous_assembly_is_not_selected_arbitrarily(tmp_path):
    opener, _ = service(ambiguous=True)
    with pytest.raises(ValueError, match="specific GCF_"):
        organisms.resolve("NC_000001", tmp_path, opener=opener)
    assert not list(tmp_path.rglob("*.json"))


@pytest.mark.parametrize(
    "failure",
    ["empty", "draft", "missing-replicon", "version", "page", "cycle", "malformed"],
)
def test_incomplete_or_mismatched_resolution_is_not_published(tmp_path, failure):
    base, _ = service()

    def opener(url):
        data = json.load(base(url))
        if failure == "empty" and "sequence_assemblies" in url:
            data = {}
        if failure == "draft" and "dataset_report" in url:
            data["reports"][0]["assembly_info"]["assembly_level"] = "Scaffold"
        if failure == "missing-replicon" and "sequence_reports" in url:
            data["reports"] = data["reports"][:1]
            data["total_count"] = 1
        if failure == "page" and "sequence_reports" in url:
            data["total_count"] = 10
        if failure == "cycle" and "sequence_reports" in url:
            data["next_page_token"] = "same"
        if failure == "malformed" and "sequence_reports" in url:
            data["reports"][0]["refseq_accession"] = "../outside"
        return response(data)

    with pytest.raises(ValueError):
        organisms.resolve(
            "NC_000001.1" if failure == "version" else "NC_000001",
            tmp_path,
            opener=opener,
        )
    assert not list(tmp_path.rglob("*.json"))


def test_failed_refresh_preserves_cache_and_successful_refresh_replaces_it(tmp_path):
    opener, _ = service()
    first = organisms.resolve("NC_000001", tmp_path, opener=opener)
    path = tmp_path / "metadata/NC_000001.json"
    before = path.read_bytes()
    with pytest.raises(ValueError):
        organisms.resolve(
            "NC_000001", tmp_path, refresh=True, opener=lambda _: response({})
        )
    assert path.read_bytes() == before

    def changed(url):
        data = json.load(opener(url))
        if "dataset_report" in url:
            data["reports"][0]["organism"]["organism_name"] = "Updated species"
        return response(data)

    second = organisms.resolve("NC_000001", tmp_path, refresh=True, opener=changed)
    assert first["organism"] != second["organism"] == "Updated species"


def test_parallel_resolvers_publish_once(tmp_path):
    opener, urls = service()
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(
            pool.map(
                lambda _: organisms.resolve("NC_000001", tmp_path, opener=opener),
                range(4),
            )
        )
    assert len(urls) == 4
    assert results.count(results[0]) == 4


def test_duplicate_organisms_through_different_replicons_are_rejected(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        organisms,
        "resolve",
        lambda query, *a, **kw: metadata(
            query, replicons=["NC_000001.2", "NC_000002.1"]
        ),
    )
    with pytest.raises(ValueError, match="Duplicate organism"):
        organisms.resolve_many(["NC_000001", "NC_000002"], tmp_path)


@pytest.mark.parametrize(
    "change",
    [
        {"query": "NC_000002"},
        {"replicons": ["../../bad"]},
        {"replicons": []},
        {"replicons": ["NC_000001.2"] * 2},
        {"assembly": "bad"},
        {"schema_version": 2},
        {"organism": ""},
    ],
)
def test_invalid_cached_metadata_is_rejected(tmp_path, change):
    folder = tmp_path / "metadata"
    folder.mkdir()
    (folder / "NC_000001.json").write_text(json.dumps(dict(metadata(), **change)))
    with pytest.raises(ValueError):
        organisms.resolve(
            "NC_000001", tmp_path, opener=lambda _: pytest.fail("network request")
        )


def test_genomes_command_prefetches_exact_replicons_and_refreshes(
    tmp_path, monkeypatch, capsys
):
    queries, downloads = [], []

    def resolve(identifier, database, **kw):
        queries.append((identifier, database, kw))
        return metadata(identifier)

    monkeypatch.setattr(organisms, "resolve", resolve)
    monkeypatch.setattr(
        ncbi,
        "genome_path",
        lambda identifier, directory, database, **kw: downloads.append(
            (identifier, database, kw)
        ),
    )
    assert (
        main(
            [
                "genomes",
                "NC_000001",
                "--genome-db",
                str(tmp_path),
                "--download",
                "--refresh",
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out) == [metadata()]
    assert queries == [("NC_000001", tmp_path, {"refresh": True})]
    assert downloads == [("NC_000001.2", tmp_path, {"refresh": True})]


def test_labels_come_from_job_metadata_or_genbank_without_catalog(tmp_path):
    from coprarna.conservation import organism_names

    (tmp_path / "genome_metadata.json").write_text(json.dumps([metadata()]))
    assert organism_names(["NC_000001", "NC_000001.2", "NC_00001"], tmp_path) == [
        "Tes_organism_strain",
        "Tes_organism_strain",
        "NC_00001",
    ]


def test_failed_genome_refresh_keeps_old_bytes(tmp_path):
    path = tmp_path / "NC_000001.2.gb.gz"
    path.write_bytes(b"previous cache")
    with pytest.raises(ValueError, match="Incomplete"):
        ncbi.download_genbank(
            "NC_000001.2", path, refresh=True, opener=lambda _: BytesIO(b"error")
        )
    assert path.read_bytes() == b"previous cache"


def test_new_accession_does_not_require_allow_list(tmp_path):
    from coprarna.pipeline import Config, _validate

    source = tmp_path / "input.fa"
    source.write_text(
        ">NC_999999.2\nACGT\n>NZ_CP999998.1\nACGT\n>GCF_123456789.1\nACGT\n"
    )
    assert _validate(source, Config(), tmp_path) == [
        "NC_999999.2",
        "NZ_CP999998.1",
        "GCF_123456789.1",
    ]


def test_legacy_cli_and_aliases_are_removed(tmp_path):
    for command in (["status"], ["helper", "--legacy", "format-fasta", "input.fa"]):
        with pytest.raises(SystemExit):
            main(command)
    assert (
        main(["helper", "--workdir", str(tmp_path), "format_fasta.pl", "input.fa"]) == 2
    )


def test_job_local_protein_codes_are_short_and_unique():
    from coprarna.homology import protein_code

    assert [protein_code(i) for i in (0, 25, 26, 26**4 - 1)] == [
        "aaaa",
        "aaaz",
        "aaba",
        "zzzz",
    ]
    assert len({protein_code(i) for i in range(1000)}) == 1000
    for index in (-1, 26**4):
        with pytest.raises(ValueError):
            protein_code(index)


def test_metadata_output_cannot_overwrite_original_fasta(tmp_path):
    from coprarna.pipeline import Config, _validate

    source = tmp_path / "genome_metadata.json"
    original = ">NC_000001\nACGT\n>NC_000002\nACGT\n>NC_000003\nACGT\n"
    source.write_text(original)
    with pytest.raises(ValueError, match="collides"):
        _validate(source, Config(), tmp_path)
    assert source.read_text() == original


def test_request_retries_transient_errors_but_not_permanent_errors(monkeypatch):
    from urllib.error import HTTPError

    calls = []

    def transient(request, timeout):
        calls.append((request, timeout))
        if len(calls) < 3:
            raise HTTPError(request, 429, "slow down", {}, None)
        return BytesIO(b"ok")

    monkeypatch.setattr(ncbi, "urlopen", transient)
    monkeypatch.setattr(ncbi.time, "sleep", lambda _: None)
    with ncbi.open_url("https://example.invalid") as result:
        assert result.read() == b"ok"
    assert len(calls) == 3 and all(timeout == 60 for _, timeout in calls)
    calls.clear()

    def permanent(request, timeout):
        calls.append(request)
        raise HTTPError(request, 404, "missing", {}, None)

    monkeypatch.setattr(ncbi, "urlopen", permanent)
    with pytest.raises(HTTPError):
        ncbi.open_url("https://example.invalid")
    assert len(calls) == 1


def test_different_cached_revisions_are_not_independent_organisms(
    tmp_path, monkeypatch
):
    records = [
        metadata(assembly="GCF_000000001.1", replicons=["NC_000001.1"]),
        metadata(assembly="GCF_000000001.2", replicons=["NC_000001.2"]),
    ]
    monkeypatch.setattr(organisms, "resolve", lambda *a, **kw: records.pop(0))
    with pytest.raises(ValueError, match="Duplicate organism"):
        organisms.resolve_many(["NC_000001.1", "NC_000001.2"], tmp_path)
