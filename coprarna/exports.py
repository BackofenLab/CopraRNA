"""Legacy enrichment exports; text formatting is intentionally preserved."""

import json
import math
from pathlib import Path
import re


def _color(groups: int, group: int, pvalue: float, solid: bool) -> str:
    hue = group * (360 / groups)
    saturation = 1 if solid or pvalue <= 1e-6 else -0.15 * math.log10(pvalue)
    hi = int(hue / 60)
    fraction = hue / 60 - hi
    low, falling, rising = (
        1 - saturation,
        1 - saturation * fraction,
        1 - saturation * (1 - fraction),
    )
    values = {
        0: (1, rising, low),
        6: (1, rising, low),
        1: (falling, 1, low),
        2: (low, 1, rising),
        3: (low, falling, 1),
        4: (rising, low, 1),
        5: (1, low, falling),
    }[hi]
    return "#" + "".join(f"{int(value * 255):02x}" for value in values)


def _json_string(value):
    from .compatibility import legacy_enabled

    return (
        '"' + value + '"' if legacy_enabled() else json.dumps(value, ensure_ascii=False)
    )


def heatmap_json(path: Path) -> str:
    """Write copraRNA.json from the retained enrichment membership table.

    The first three lines define clusters, dimensions and feature labels. Each
    remaining membership becomes a colored JSON cell; corrected production escapes
    labels safely while an explicit historical helper context preserves old bytes.
    """
    lines = path.read_text().splitlines()
    groups, parameters, header = (line.split("\t") for line in lines[:3])
    group_text = ";".join(
        f"{i + 1}:{float(value):.2f}:{_color(int(parameters[1]), i + 1, 0, True)}"
        for i, value in enumerate(groups)
    )
    output = []
    for row_number, line in enumerate(lines[3:], 1):
        row = line.split("\t")
        for feature_index, i in enumerate(range(11, len(row)), 1):
            color = (
                _color(int(parameters[1]), int(parameters[i - 9]), float(row[1]), False)
                if float(row[i]) == 1
                else "#ffffff"
            )
            # The original emits hand-built JSON, including its escaping behavior.
            output.append(
                "{"
                + f'"flagF":{int(row_number == 1)},"fMax":{len(parameters) - 2},'
                + f'"fNum":{feature_index},"feature":{_json_string(header[i])},'
                + f'"flagS":{int(feature_index == 1)},"sMax":{parameters[0]},'
                + f'"sNum":{row_number},"source":{_json_string(row[4])},'
                + f'"color":"{color}","value":{row[1]},"group":"{group_text}"'
                + "}"
            )
    (path.parent / "copraRNA.json").write_text(
        "[\n" + ",\n".join(output) + ("\n" if output else "") + "]\n"
    )
    return ""


def specific_targets(directory: Path) -> str:
    """Return enriched-term targets found only in the organism-specific prediction.

    Compare DAVID cluster/chart reports with the annotated IntaRNA table, retaining
    the original GeneID lookup and output field ordering for helper parity.
    """
    main, single = {}, {}
    score = 0.0
    for line in (directory / "termClusterReport.txt").read_text().splitlines():
        match = re.search(r"Enrichment\s+Score:\s+(\d+\.\d+)", line)
        if match:
            score = float(match[1])
        fields = line.split("\t")
        if len(fields) > 5 and fields[1] != "Term" and score >= 1:
            main[fields[1]] = fields[5]
    lines = (directory / "IntaRNA_chartReport_grepped.txt").read_text().splitlines()
    for i in range(0, len(lines), 2):
        single[lines[i + 1].split("=")[1]] = lines[i].split("=")[1]
    intarna = (directory / "intarna_websrv_table_ncbi.csv").read_text().splitlines()
    output = [
        "term;predicted_targets EntrezGeneID(locus_tag|gene_name|target_start|target_end)\n"
    ]
    for term in main:
        if term not in single:
            continue
        known = {re.sub(r"\s+", "", value) for value in main[term].split(",")}
        extra = []
        for value in single[term].split(","):
            gid = re.sub(r"\s+", "", value)
            if gid not in known:
                # Legacy grep performs substring/regex matching, including multiple hits.
                joined = "\n".join(line for line in intarna if re.search(gid, line))
                fields = joined.split(";")
                extra.append(f"{gid}({fields[0]}|{fields[36]}|{fields[8]}|{fields[9]})")
        if extra:
            output.append(term + ";" + ";".join(extra) + "\n")
    return "".join(output)
