"""Service request parity with offline SOAP/HTTP doubles and native rendering."""

from contextlib import redirect_stdout
from datetime import datetime, timezone
from io import BytesIO, StringIO
import os
import runpy
import sys
from types import SimpleNamespace
from urllib.parse import urlparse, parse_qs

import pytest

from coprarna import david, maintenance, ncbi, webserver
from tests.support.legacy import ROOT, FIXTURES


@pytest.mark.parametrize("copra", [True, False])
def test_david_requests_and_stdout_match_untouched_source(copra, tmp_path, monkeypatch):
    pytest.importorskip("suds.client")
    import suds.client

    calls = []

    class Service:
        def __getattr__(self, name):
            def method(*args):
                calls.append((name, args))
                return name + " result"

            return method

    class Client:
        def __init__(self, url):
            calls.append(("Client", (url,)))
            self.service = Service()
            self.wsdl = SimpleNamespace(
                services=[
                    SimpleNamespace(
                        setlocation=lambda url: calls.append(("setlocation", (url,)))
                    )
                ]
            )

        def __str__(self):
            return "offline DAVID client"

    monkeypatch.setattr(suds.client, "Client", Client)
    path = tmp_path / "prediction.csv"
    if copra:
        path.write_bytes(
            (FIXTURES / "functional-enrichment/input/prediction.csv").read_bytes()
        )
    else:
        path.write_text(
            ";".join(["header"] * 39)
            + "\n"
            + "".join(
                ";".join(["x"] * 37 + [str(i), "annotation"]) + "\n" for i in range(10)
            )
        )
    source = (
        "DAVIDWebService_CopraRNA.py"
        if copra
        else "DAVIDWebService_IntaRNA_chartReport.py"
    )
    monkeypatch.setattr(sys, "argv", [source, str(path), "5"])
    monkeypatch.setattr(sys, "path", list(sys.path))
    out = StringIO()
    with redirect_stdout(out):
        runpy.run_path(str(ROOT / "coprarna_aux" / source), run_name="__main__")
    expected_calls = list(calls)
    calls.clear()
    actual = (david.coprarna_report if copra else david.intarna_report)(path, 5, Client)
    assert actual == out.getvalue()
    assert calls == expected_calls


def test_ncbi_requests_first_hit_and_preserves_download_bytes(tmp_path):
    responses = [
        b"<eSearchResult><IdList><Id>123</Id><Id>456</Id></IdList></eSearchResult>",
        b"<eSummaryResult/>",
        b"LOCUS cached\n//\n",
    ]
    urls = []

    def opener(url):
        urls.append(url)
        return BytesIO(responses.pop(0))

    path = tmp_path / "genome.gb"
    ncbi.download_genbank("NC_000913", path, opener=opener)
    assert path.read_bytes() == b"LOCUS cached\n//\n"
    requests = [
        (urlparse(url).path.rsplit("/", 1)[-1], parse_qs(urlparse(url).query))
        for url in urls
    ]
    assert [name for name, _ in requests] == [
        "esearch.fcgi",
        "esummary.fcgi",
        "efetch.fcgi",
    ]
    assert requests[0][1]["term"] == ["NC_000913"]
    assert requests[2][1] == dict(
        db=["nuccore"],
        id=["123"],
        tool=["intarna_web"],
        email=["intarna@informatik.uni-freiburg.de"],
        rettype=["gbwithparts"],
    )
    ncbi.download_genbank(
        "NC_000913", path, opener=lambda _: pytest.fail("cached file downloaded again")
    )


def test_missing_ncbi_accession(tmp_path):
    with pytest.raises(ValueError, match="no genome"):
        ncbi.download_genbank(
            "absent",
            tmp_path / "none.gb",
            opener=lambda _: BytesIO(b"<eSearchResult><IdList/></eSearchResult>"),
        )


def test_maintenance_mapping_matches_seeded_legacy(tmp_path):
    fixture = FIXTURES / "maintenance"
    maintenance.build_organism_mapping(
        tmp_path,
        fixture / "prokaryotes.txt",
        seed=123,
        now=datetime(2020, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
    )
    for name in ("kegg2refseqnew.csv", "CopraRNA_available_organisms.txt"):
        assert (tmp_path / name).read_bytes() == (fixture / name).read_bytes()
    assert not (tmp_path / "CopraRNA_available_organisms.tmp").exists()


@pytest.mark.skipif(
    os.environ.get("COPRARNA_BROWSER_TESTS") != "1", reason="opt-in Chromium rendering"
)
def test_rasterize_real_chromium(tmp_path):
    html = tmp_path / "page.html"
    html.write_text(
        '<html><body><h1>CopraRNA report</h1><svg width="50" height="50"><rect width="50" height="50" fill="green"/></svg></body></html>'
    )
    for extension in ("png", "pdf"):
        output = tmp_path / ("report." + extension)
        webserver.rasterize_heatmap(str(html), output)
        assert output.stat().st_size > 500
        assert output.read_bytes().startswith(
            b"%PDF" if extension == "pdf" else b"\x89PNG"
        )
