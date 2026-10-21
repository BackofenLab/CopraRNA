"""NCBI accession lookup and GenBank download without BioPerl."""

from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def download_genbank(accession: str, path: Path, *, opener=None) -> str:
    if path.exists():
        return ""
    opener = opener or urlopen

    def request(endpoint, parameters):
        with opener(BASE + endpoint + "?" + urlencode(parameters)) as response:
            return response.read()

    search = request(
        "esearch.fcgi",
        dict(db="nuccore", field="Accession", term=accession, tool="coprarna"),
    )
    identifier = ET.fromstring(search).findtext("./IdList/Id")
    if not identifier:
        raise ValueError(f"no genome with accession number '{accession}' found")
    request("esummary.fcgi", dict(db="nuccore", id=identifier))
    data = request(
        "efetch.fcgi",
        dict(
            db="nuccore",
            id=identifier,
            tool="intarna_web",
            email="intarna@informatik.uni-freiburg.de",
            rettype="gbwithparts",
        ),
    )
    path.write_bytes(data)
    return ""
