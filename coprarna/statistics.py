"""Validated deterministic statistical preprocessing from the legacy R scripts.

GEV fitting and correlated p-value joining have separate migration requirements;
no substitute numerical algorithm is used here.
"""

from decimal import Decimal
import math
from pathlib import Path


def bh_adjust(pvalues: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg correction to finite positions, preserving missing NA."""
    valid = [i for i, value in enumerate(pvalues) if not math.isnan(value)]
    ordered = sorted(valid, key=lambda i: pvalues[i], reverse=True)
    adjusted = [math.nan] * len(pvalues)
    current = 1.0
    for rank, index in enumerate(ordered):
        current = min(current, len(valid) / (len(valid) - rank) * pvalues[index])
        adjusted[index] = current
    return adjusted


def r_number(value: float) -> str:
    """Format a numeric scalar with R's retained precision and NA/Inf conventions."""
    if math.isnan(value):
        return "NA"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    if value == 0:
        return "0"
    fixed = format(Decimal(format(value, ".15g")), "f")
    if "." in fixed:
        fixed = fixed.rstrip("0").rstrip(".")
    mantissa, exponent = format(value, ".14e").split("e")
    scientific = mantissa.rstrip("0").rstrip(".") + "e" + exponent
    return scientific if len(scientific) < len(fixed) else fixed


def adjust_pvalues(path: Path, directory: Path) -> str:
    """Read the headered p-value table and write the BH-adjusted padj.csv column."""
    rows = [line.split() for line in path.read_text().splitlines()[1:] if line.strip()]
    width = len(rows[0]) if rows else 0
    values = [
        float(row[col]) if row[col] != "NA" else math.nan
        for col in range(width)
        for row in rows
    ]
    adjusted = bh_adjust(values)
    (directory / "padj.csv").write_text(
        "fdr\n" + "".join(r_number(p) + "\n" for p in adjusted)
    )
    return ""


def _r_normalize_numeric_columns(rows: list[list[str]]) -> list[list[str]]:
    # read.table type-converts each column before write.table prints it.
    """Reproduce read.table type conversion followed by write.table formatting."""
    if not rows:
        return rows
    result = [list(row) for row in rows]
    for col in range(len(rows[0])):
        try:
            values = [
                float(row[col]) if row[col] not in ("NA", "") else math.nan
                for row in rows
            ]
        except ValueError:
            continue
        for row, value in zip(result, values):
            row[col] = r_number(value)
    return result


def prepare_copra2(path: Path, directory: Path) -> str:
    """Write ordered homolog tag sets for the comparative cluster-joining stage."""
    rows = _r_normalize_numeric_columns(
        [line.split(";") for line in path.read_text().splitlines()[1:]]
    )
    groups = dict.fromkeys(row[-1] for row in rows)
    output = []
    for cluster in groups:
        members = sorted(
            (row for row in rows if row[-1] == cluster), key=lambda row: row[1]
        )
        output.append("0;" + ";".join(row[0] for row in members) + ";\n")
    (directory / "CopraRNA2_prep.csv").write_text("".join(output))
    return ""
