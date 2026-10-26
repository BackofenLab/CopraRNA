"""Literal ports of the text and IntaRNA table helpers."""

from collections import defaultdict
import math
from pathlib import Path
import re
from .io import paths, read_text, write_text

INTARNA_HEADER = "id1;id2;seq1;seq2;subseq1;subseq2;subseqDP;subseqDB;start1;end1;start2;end2;hybridDP;hybridDB;E;ED1;ED2;Pu1;Pu2;E_init;E_loops;E_dangleL;E_dangleR;E_endL;E_endR;seedStart1;seedEnd1;seedStart2;seedEnd2;seedE;seedED1;seedED2;seedPu1;seedPu2;E_norm"


def perl_split(pattern: str, value: str) -> list[str]:
    """Split with a regular expression and discard trailing empty fields as Perl does."""
    parts = re.split(pattern, value)
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def blast_to_homology(path: Path) -> str:
    """Return DomClust similarities from BLAST tabular format 6, plain or gzip.

    Canonicalize each protein pair, keep its best eligible e-value, and reproduce
    the original identity-to-distance transform and score scaling. This is the
    interface adapter shared by historical fixtures and the current BLAST+ caller.
    """
    rows = []
    for line in read_text(path).splitlines():
        if line.startswith("#") or not line.strip():
            continue
        row = line.split()
        if row[0] > row[1]:
            row = [row[1], row[0], *row[2:6], *row[8:10], *row[6:8], *row[10:]]
            line = "\t".join(row)
        rows.append((row, line))
    rows.sort(key=lambda pair: (pair[0][0], pair[0][1], float(pair[0][10]), pair[1]))
    previous = None
    output = []
    for row, _ in rows:
        if float(row[10]) > 0.001 or tuple(row[:2]) == previous:
            continue
        previous = tuple(row[:2])
        distance = 100 - float(row[2])
        fraction = distance / 100
        term = 1 - fraction - 0.2 * fraction * fraction
        if term <= 0:
            continue
        if distance > 0:
            distance = -math.log(term) * 100
        score = float(row[11]) / (1 / 3)
        output.append(
            " ".join([*row[:2], *row[6:10], f"{distance:.0f}", f"{score:.0f}"]) + "\n"
        )
    return "".join(output)


def disentangle(directory: Path) -> str:
    """Write optimal/suboptimal IntaRNA tables from sorted native predictions.

    Keep the best and second-best interaction per target. Read/write gzip
    directly when the input is compressed and keep the required native schema.
    """
    for path in paths(directory, "*.fa.intarna.sorted.csv"):
        optimal, suboptimal = {}, {}
        for line in read_text(path).splitlines()[1:]:
            tag = line.split(";")[0]
            if tag in optimal:
                suboptimal.setdefault(tag, line)
            else:
                optimal[tag] = line
        for suffix, rows in [("opt", optimal), ("subopt", suboptimal)]:
            name = (
                path.name.split(".")[0]
                + f"_{suffix}.intarna.csv"
                + (".gz" if path.suffix == ".gz" else "")
            )
            write_text(
                directory / name,
                INTARNA_HEADER + "\n" + "".join(line + "\n" for line in rows.values()),
            )
    return ""


def cluster_predictions(directory: Path) -> str:
    """Return the best eligible organism interaction for each protein-homology cluster.

    Read cluster.tab and optimal p-value tables, remove DomClust domain suffixes
    for locus lookup, and append the original cluster number to each selected row.
    """
    lookup = {}
    for path in paths(directory, "*_opt.intarna.csv"):
        for line in read_text(path).splitlines()[1:]:
            line = line.lower()
            lookup[line.split(";")[0]] = line
    output = [INTARNA_HEADER[1:] + ";p-value;clusternumber\n"]
    for index, line in enumerate(
        (directory / "cluster.tab").read_text().splitlines()[1:], 1
    ):
        for cell in perl_split(r"\t", line.lower())[1:]:
            pvalue, chosen = 100.0, ""
            for entry in re.split(r"\s", cell):
                tag = re.sub(r"\w+:", "", re.sub(r"\(\d+\)", "", entry))
                if tag in lookup:
                    fields = lookup[tag].split(";")
                    candidate = float(fields[35])
                    if candidate < pvalue:
                        chosen, pvalue = lookup[tag], candidate
            if chosen:
                output.append(f"{chosen};{index}\n")
    return "".join(output)


def add_sampled_count(path: Path) -> str:
    """Append the reserved Amount sampled report column, always zero for CopraRNA2.

    Keep the table/frontend schema without retaining obsolete missing-value sampling.
    """
    lines = read_text(path).splitlines()
    output = [lines[0] + ",Amount sampled\n"]
    for line in lines[1:]:
        count = 0
        output.append(f"{line},{count}\n")
    return "".join(output)


def additional_homologs(path: Path, clusters: Path) -> str:
    """Append same-organism homolog tags from cluster.tab to each annotated row.

    Preserve first-seen order, suppress duplicates and omit the selected target.
    """
    lookup = defaultdict(list)
    for line in clusters.read_text().splitlines():
        if line.startswith(("#", "running")):
            continue
        for cell in perl_split(r"\t", line)[1:]:
            homologs = []
            for entry in re.split(r"\s", cell):
                match = re.search(r"^[^:]+:(.*)", entry)
                if match:
                    tag = re.sub(r"\(\d+\).*", "", match[1]).lower()
                    homologs.append(tag)
            for tag in homologs:
                lookup[tag].extend(homologs)
    lines = read_text(path).splitlines()
    output = [lines[0] + ",Additional homologs\n"]
    for line in lines[1:]:
        fields = line.split(",")
        match = (
            re.match(r"^(\S+?)\(\S+GeneID\S+\)", fields[1]) if len(fields) > 1 else None
        )
        extra = ""
        if match:
            extra = "".join(
                tag + " " for tag in dict.fromkeys(lookup[match[1]]) if tag != match[1]
            )
        output.append(line + "," + extra + "\n")
    return "".join(output)


def enrichment_report(path: Path) -> str:
    """Return a DAVID annotation-cluster table from the SOAP report text.

    Parse named fields and scores, retain the original column order and strip
    embedded quotes to match the committed service-adapter oracle.
    """
    lines = read_text(path).splitlines()
    header = "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n"
    fields = [
        "categoryName",
        "termName",
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
    output = []
    count = 0
    for i, line in enumerate(lines):
        if "score=" in line:
            score = line.split("=")[1]
            output.append(
                ("\n" if count else "")
                + f"Annotation Cluster {count + 1}\tEnrichment Score:  {score}\n"
                + header
            )
            count += 1
        if "afdr=" in line:
            values = {}
            for j, entry in enumerate(lines[i:]):
                if j > 0 and "afdr=" in entry:
                    break
                pair = entry.split("=")
                values[pair[0]] = pair[1] if len(pair) > 1 else ""
            output.append(
                "\t".join(values.get(key, "") for key in fields).replace('"', "") + "\n"
            )
    return "".join(output)
