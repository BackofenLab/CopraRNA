"""Python 3 DAVID SOAP adapters retaining legacy requests and report formatting."""

from io import StringIO
from pathlib import Path

WSDL = "https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl"
ENDPOINT = "https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/"


def _report(path, count, copra, client_factory=None):
    lines = path.read_text().splitlines(keepends=True)
    background = []
    for line in lines[1:]:
        fields = line.split("," if copra else ";")
        if copra:
            if fields[2]:
                background.append(fields[2].split(":")[1][:-1])
        else:
            background.append(fields[37])
    if client_factory is None:
        from suds.client import Client

        client_factory = Client
    output = StringIO()

    def emit(value):
        print(value, file=output)

    emit("background:" + str(len(background)))
    emit("input:" + str(len(background[:count])))
    emit("url=" + WSDL)
    client = client_factory(WSDL)
    client.wsdl.services[0].setlocation(ENDPOINT)
    emit(client)
    emit(client.service.authenticate("patrickrw@gmx.net"))
    emit(
        client.service.addList(
            ",".join(background[:count]), "ENTREZ_GENE_ID", "make_up", 0
        )
    )
    if copra:
        emit(client.service.addList(",".join(background), "ENTREZ_GENE_ID", "back", 1))
    emit(client.service.getDefaultCategoryNames())
    emit(
        client.service.getTermClusterReport(3, 2, 2, 0.5, 85)
        if copra
        else client.service.getChartReport(1, 1)
    )
    return output.getvalue()


def coprarna_report(path: Path, count: int, client_factory=None) -> str:
    """Request DAVID term clusters for the first count comparative predictions.

    Use all supplied predictions as background and return the captured request/
    response transcript. An injected SOAP client enables offline request parity.
    """
    return _report(path, count, True, client_factory)


def intarna_report(path: Path, count: int, client_factory=None) -> str:
    """Request a DAVID chart for the first count annotated native interactions.

    Return the service transcript; an injected client replaces network access in
    tests without changing endpoint, account or category/request conventions.
    """
    return _report(path, count, False, client_factory)
