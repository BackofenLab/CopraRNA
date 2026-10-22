"""Service request parity with offline SOAP/HTTP doubles and native rendering."""

from io import BytesIO
import os
import json
from types import SimpleNamespace
from urllib.parse import urlparse, parse_qs

import pytest

from coprarna import david, ncbi, webserver
from tests.support.fixtures import FIXTURES


@pytest.mark.parametrize("copra", [True, False])
def test_david_requests_and_stdout_match_captured_source(copra, tmp_path, monkeypatch):
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
    expected = json.loads(
        (FIXTURES / "david" / ("copra.json" if copra else "intarna.json")).read_text()
    )
    actual = (david.coprarna_report if copra else david.intarna_report)(path, 5, Client)
    assert actual == expected["stdout"]
    assert json.loads(json.dumps(calls)) == expected["calls"]


def test_ncbi_requests_exact_accession_and_preserves_download_bytes(tmp_path):
    responses = [
        b"LOCUS cached\n//\n",
    ]
    urls = []

    def opener(url):
        urls.append(getattr(url, "full_url", url))
        return BytesIO(responses.pop(0))

    path = tmp_path / "genome.gb.gz"
    ncbi.download_genbank("NC_000913", path, opener=opener)
    import gzip

    assert gzip.decompress(path.read_bytes()) == b"LOCUS cached\n//\n"
    requests = [
        (urlparse(url).path.rsplit("/", 1)[-1], parse_qs(urlparse(url).query))
        for url in urls
    ]
    assert requests == [
        (
            "efetch.fcgi",
            dict(
                db=["nuccore"],
                id=["NC_000913"],
                tool=["coprarna"],
                rettype=["gbwithparts"],
                retmode=["text"],
            ),
        )
    ]
    ncbi.download_genbank(
        "NC_000913", path, opener=lambda _: pytest.fail("cached file downloaded again")
    )


def test_missing_ncbi_accession(tmp_path):
    with pytest.raises(ValueError, match="Incomplete"):
        ncbi.download_genbank(
            "absent",
            tmp_path / "none.gb.gz",
            opener=lambda _: BytesIO(b"<eSearchResult><IdList/></eSearchResult>"),
        )


@pytest.mark.skipif(
    os.environ.get("COPRARNA_BROWSER_TESTS") != "1", reason="opt-in Chromium rendering"
)
def test_rasterize_real_chromium(tmp_path):
    html = tmp_path / "page.html"
    html.write_text(
        '<html><body><h1>CopraRNA report</h1><svg width="50" height="50"><rect width="50" height="50" fill="green"/></svg></body></html>'
    )
    for extension in ("png", "pdf", "svg"):
        output = tmp_path / ("report." + extension)
        webserver.rasterize_heatmap(str(html), output)
        if extension == "svg":
            import xml.etree.ElementTree as ET

            drawing = ET.parse(output).getroot()
            assert any(node.tag.endswith("rect") for node in drawing.iter())
            assert not any(node.tag.endswith("image") for node in drawing.iter())
        else:
            assert output.stat().st_size > 500
        assert output.read_bytes().startswith(
            b"%PDF"
            if extension == "pdf"
            else b"<?xml"
            if extension == "svg"
            else b"\x89PNG"
        )
