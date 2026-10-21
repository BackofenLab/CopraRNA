"""DAVID enrichment processing and the original heatmap table format."""

from pathlib import Path
import re
from importlib.resources import files
from .rcompat import read_table, write_table
from .statistics import r_number
from . import david, tables, exports, annotation


def extract_enriched_terms(predictions: Path, report: Path, output: Path) -> str:
    rows = [line.split("\t") for line in report.read_text().splitlines() if line]
    starts = [
        i for i, row in enumerate(rows) if len(row) > 1 and "Enrichment Score" in row[1]
    ]
    scores = [float(rows[i][1].replace("Enrichment Score: ", "")) for i in starts]
    count = max(1, sum(score >= 1 for score in scores))
    groups = []
    for i in range(min(count, len(starts))):
        group = rows[
            starts[i] + 2 : starts[i + 1] if i + 1 < len(starts) else len(rows)
        ]
        seen = set()
        filtered = []
        for row in group:
            if row[0] == "KEGG_PATHWAY":
                pathway = row[1].split(":")[1]
                if pathway in seen:
                    continue
                seen.add(pathway)
            filtered.append(row)
        groups.append(filtered)
    identifiers = list(
        dict.fromkeys(
            int(gid.strip())
            for group in groups
            for row in group
            for gid in row[5].split(",")
        )
    )
    header, prediction_rows = read_table(predictions)
    terms = [
        r_number(round(float(row[9]), 2)) + "  " + row[1]
        for group in groups
        for row in group
    ]
    memberships = [j + 1 for j, group in enumerate(groups) for row in group]
    output_header = [
        "Entrez GeneID",
        "CopraRNA_pValue",
        "Locustag",
        "GeneName",
        "",
        "IntaRNA energy",
        "IntaRNA p-Value",
        "start mRNA",
        "end mRNA",
        "start ncRNA",
        "end ncRNA",
        *terms,
    ]
    body = []
    for identifier in identifiers:
        match = next(
            (row for row in prediction_rows if re.search(str(identifier), row[2])), None
        )
        if match is None:
            continue
        parts = match[2].split("|")
        name = parts[0].split("(")
        locus = name[0] if len(name) == 2 else "NA"
        gene = name[1] if len(name) == 2 else "NA"
        label = locus + " - " + gene
        presence = [
            int(bool(re.search(str(identifier), row[5])))
            for group in groups
            for row in group
        ]
        body.append(
            [
                str(identifier),
                r_number(float(match[1])),
                locus,
                label,
                label,
                *parts[1:7],
                *map(str, presence),
            ]
        )
    width = len(output_header)
    meta2 = [r_number(score) for score in scores[: len(groups)]] + [""] * (
        width - len(groups)
    )
    meta = [str(len(body)), str(len(groups)), *map(str, memberships), *([""] * 9)]
    write_table(output, None, [meta2, meta, output_header, *body], sep="\t")
    return ""


def _grep_report(text, fields):
    output = []
    for line in text.splitlines():
        if any(re.search(field + r"\s=", line) for field in fields):
            output.append(line.lstrip().replace(" = ", "=").replace(", ", ","))
    return "\n".join(output) + "\n"


def enrich_prediction(directory, count):
    from .webserver import rasterize_heatmap

    try:
        report = david.coprarna_report(directory / "CopraRNA_result_all.csv", count)
    except Exception:
        report = ""  # The legacy controller continues after a failed SOAP subprocess.
    (directory / "DAVID_enrichment_temp.txt").write_text(report)
    fields = [
        "termName",
        "categoryName",
        "score",
        "listHits",
        "percent",
        "ease",
        "geneIds",
        "listTotals",
        "popHits",
        "popTotals",
        "foldEnrichment",
        "bonferroni",
        "benjamini",
        "afdr",
    ]
    grepped = _grep_report(report, fields)
    (directory / "DAVID_enrichment_grepped_temp.txt").write_text(grepped)
    term = tables.enrichment_report(directory / "DAVID_enrichment_grepped_temp.txt")
    if not term:
        term = "If you are reading this, then your prediction did not return an enrichment, your organism of interest is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your CopraRNA\nprediction or create your enrichment manually at the DAVID homepage.\n"
    (directory / "termClusterReport.txt").write_text(term)
    if "Annotation Cluster" in term:
        extract_enriched_terms(
            directory / "CopraRNA_result_all.csv",
            directory / "termClusterReport.txt",
            directory / "enrichment.txt",
        )
        exports.heatmap_json(directory / "enrichment.txt")
        for name in ("copra_heatmap.html", "index-thumb.html", "index-pdf.html"):
            (directory / name).write_bytes(
                files("coprarna").joinpath("data/" + name).read_bytes()
            )
        rasterize_heatmap(
            str(directory / "index-thumb.html"), directory / "enriched_heatmap_big.png"
        )
        rasterize_heatmap(
            str(directory / "index-pdf.html"), directory / "enriched_heatmap_big.pdf"
        )
        for name in ("index-thumb.html", "index-pdf.html"):
            (directory / name).unlink()
    # The auxiliary single-organism report uses the same requested top count.
    from .runtime import options
    from .intarna import numeric_sort
    from Bio import SeqIO

    opt = options(directory)
    ooi = next(SeqIO.parse(directory / "ncrna.fa", "fasta")).id.removeprefix("ncRNA_")
    path = (
        directory
        / f"{ooi}_upfromstartpos_{opt['nt upstream']}_down_{opt['nt downstream']}_opt.intarna.csv"
    )
    (directory / "intarna_websrv_table.csv").write_text(
        "".join(numeric_sort(path.read_text().splitlines(keepends=True), 35))
    )
    annotation.annotate_intarna(directory)
    try:
        chart = david.intarna_report(directory / "intarna_websrv_table_ncbi.csv", count)
    except Exception:
        chart = ""
    (directory / "IntaRNA_chartReport.txt").write_text(chart)
    (directory / "IntaRNA_chartReport_grepped.txt").write_text(
        _grep_report(chart, ["geneIds", "termName"]).replace('"', "")
    )
    specific = exports.specific_targets(directory)
    (directory / "org_of_interest_aux_enrichment.txt").write_text(specific)
    _, intarna = read_table(directory / "intarna_websrv_table_ncbi.csv", sep=";")
    terms = specific.splitlines()[1:]
    gids = sorted(set(re.findall(r";(\d+)\(", specific)))
    output = [
        "locus_tag,start_tar,stop_tar,start_query,stop_query,energy,p-value,gene_name,gene_id,annotation,functional_terms\n"
    ]
    for gid in gids:
        matches = [row for row in intarna if row[37] == gid]
        for row in matches:
            values = [row[i] for i in (0, 8, 9, 10, 11, 14, 35, 36, 37, 38)]
            labels = ";".join(line.split(";")[0] for line in terms if gid in line)
            output.append(",".join(values) + "," + labels.replace(",", " ", 1) + "\n")
    (directory / "aux_table.csv").write_text("".join(output))
