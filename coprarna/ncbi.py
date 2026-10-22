"""Streaming NCBI downloads and a reusable compressed GenBank database."""

from collections import deque
from contextlib import contextmanager
import gzip
from pathlib import Path
import shutil
import tempfile
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def download_genbank(accession: str, path: Path, *, opener=None) -> str:
    """Stream one NCBI GenBank accession into an atomically published gzip file.

    Reuse an existing path. Request HTTP gzip, preserve it if supplied, otherwise
    compress during download. Validate the terminating record marker and gzip
    integrity by streaming; incomplete responses never become cache entries.
    ``opener`` permits deterministic request tests without a live service."""
    path = Path(path)
    if path.exists():
        return ""
    if path.suffix != ".gz":
        raise ValueError("GenBank downloads require a .gz output path")
    opener = opener or urlopen

    def url(endpoint, **values):
        return BASE + endpoint + "?" + urlencode(values)

    with opener(
        url(
            "esearch.fcgi",
            db="nuccore",
            field="Accession",
            term=accession,
            tool="coprarna",
        )
    ) as response:
        identifier = ET.fromstring(response.read()).findtext("./IdList/Id")
    if not identifier:
        raise ValueError(f"no genome with accession number '{accession}' found")
    with opener(url("esummary.fcgi", db="nuccore", id=identifier)) as response:
        response.read()
    request = Request(
        url(
            "efetch.fcgi",
            db="nuccore",
            id=identifier,
            tool="intarna_web",
            email="intarna@informatik.uni-freiburg.de",
            rettype="gbwithparts",
        ),
        headers={"Accept-Encoding": "gzip"},
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix="." + path.name, dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
    try:
        with opener(request) as response:
            compressed = (
                getattr(response, "headers", {}).get("Content-Encoding", "").lower()
                == "gzip"
            )
            if compressed and path.suffix == ".gz":
                with temporary.open("wb") as target:
                    shutil.copyfileobj(response, target, 1024 * 1024)
            else:
                source = gzip.GzipFile(fileobj=response) if compressed else response
                output = (
                    gzip.open(temporary, "wb")
                    if path.suffix == ".gz"
                    else temporary.open("wb")
                )
                with output as target:
                    shutil.copyfileobj(source, target, 1024 * 1024)
        # Truncated downloads never become visible as cached genomes.
        with (
            gzip.open(temporary, "rt") if path.suffix == ".gz" else temporary.open()
        ) as stream:
            tail = deque((line.strip() for line in stream if line.strip()), maxlen=1)
        if not tail or tail[-1] != "//":
            raise ValueError(f"Incomplete GenBank download for {accession}")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return ""


@contextmanager
def _lock(path):
    """Hold a Linux advisory lock for one accession across threads and processes."""
    import fcntl

    with path.open("a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def genome_path(accession, directory, database=None):
    """Resolve a compressed genome without making a per-job copy.

    Search the shared database, then supplied job files. Preserve plain supplied
    files and stream a gzip cache copy once. Serialize cache misses per accession
    so concurrent jobs share one completed download and never read partial data."""
    directory = Path(directory)
    database = (
        Path(database).expanduser().resolve() if database else directory / "genomes"
    )
    database.mkdir(parents=True, exist_ok=True)
    target = database / (accession + ".gb.gz")
    with _lock(database / (accession + ".lock")):
        for suffix in (".gb.gz", ".gbff.gz", ".gbk.gz"):
            cached = database / (accession + suffix)
            if cached.is_file():
                return cached
        for suffix in (".gb.gz", ".gbff.gz", ".gbk.gz"):
            supplied = directory / (accession + suffix)
            if supplied.is_file():
                return supplied.resolve()
        plain = next(
            (
                directory / (accession + suffix)
                for suffix in (".gb", ".gbff", ".gbk")
                if (directory / (accession + suffix)).is_file()
            ),
            None,
        )
        if plain is not None:
            with tempfile.NamedTemporaryFile(
                prefix="." + target.name, dir=database, delete=False
            ) as stream:
                temporary = Path(stream.name)
            try:
                with plain.open("rb") as source, gzip.open(temporary, "wb") as output:
                    shutil.copyfileobj(source, output, 1024 * 1024)
                temporary.replace(target)
            finally:
                temporary.unlink(missing_ok=True)
        else:
            download_genbank(accession, target)
    return target
