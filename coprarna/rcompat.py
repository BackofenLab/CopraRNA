"""R table serialization and persisted matrix conventions used by CopraRNA."""

import csv
import math
import re
import numpy as np
from .statistics import r_number
from .io import open_text, write_text


def read_table(path, sep=","):
    """Read a header and data rows from a plain/gzip delimited table."""
    with open_text(path) as stream:
        rows = list(csv.reader(stream, delimiter=sep))
    return rows[0], rows[1:]


def make_names(names):
    """Normalize report field punctuation as required by retained R column access."""
    return [re.sub(r"[^\w.]", ".", name) for name in names]


def scalar(value, quote=False):
    """Serialize one R-compatible table value, including booleans and missing data."""
    if value is None or isinstance(value, (float, np.floating)) and math.isnan(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (int, np.integer)):
        return str(value)
    if isinstance(value, (float, np.floating)):
        return r_number(float(value))
    value = str(value)
    return '"' + value.replace('"', '\\"') + '"' if quote else value


def write_table(path, header, rows, sep=",", quote=False, row_names=False):
    """Write retained R table formatting, with optional row names and gzip storage."""
    lines = []
    if header is not None:
        lines.append(sep.join(scalar(v, quote) for v in header))
    for i, row in enumerate(rows, 1):
        cells = (
            [str(i) if row_names is True else str(row_names[i - 1])]
            if row_names
            else []
        ) + list(row)
        lines.append(sep.join(scalar(v, quote) for v in cells))
    write_text(path, "\n".join(lines) + "\n")


def numeric_columns(rows):
    """Convert wholly numeric table columns, representing empty/NA cells as NaN."""
    result = [list(row) for row in rows]
    if not rows:
        return result
    for j in range(len(rows[0])):
        try:
            values = [
                float(row[j]) if row[j] not in ("", "NA") else math.nan for row in rows
            ]
        except (TypeError, ValueError):
            continue
        for row, value in zip(result, values):
            row[j] = value
    return result


class RMatrix:
    """Carry matrix values and dimnames for R-compatible column-major serialization."""

    def __init__(self, values, columns=None, rows=None):
        self.values = np.asarray(values, dtype=object)
        self.columns = columns
        self.rows = rows


class RNamedList:
    """Carry list members and names for R-compatible serialization."""

    def __init__(self, values, names):
        self.values = values
        self.names = names


def save_rdata(path, name, values, *, names=None, columns=None):
    """Persist named lists/matrices with R dimensions, names and column-major order."""
    import rdata

    def constructor(matrix, converter):
        data = converter.convert_to_r_object(matrix.values.ravel(order="F"))
        attributes = {"dim": np.array(matrix.values.shape, dtype=np.int32)}
        if matrix.columns is not None or matrix.rows is not None:
            attributes["dimnames"] = [
                None if matrix.rows is None else np.array(matrix.rows, dtype="U"),
                None if matrix.columns is None else np.array(matrix.columns, dtype="U"),
            ]
        data.attributes = converter.convert_to_r_attributes(attributes)
        data.info.attributes = True
        return data

    def named_list(value, converter):
        data = converter.convert_to_r_object(value.values)
        data.attributes = converter.convert_to_r_attributes(
            {"names": np.array(value.names, dtype=object)}
        )
        data.info.attributes = True
        return data

    if columns is not None:
        values = [RMatrix(value, columns) for value in values]
    rdata.write_rda(
        path,
        {name: values},
        constructor_dict={RMatrix: constructor, RNamedList: named_list},
    )


def load_rdata(path, name):
    """Load the requested R object while restoring character-matrix dimensions."""
    import rdata
    from rdata.parser import RObjectType

    parsed = rdata.parser.parse_file(path)

    def convert(obj):
        if obj.info.type == RObjectType.VEC:
            return [convert(value) for value in obj.value]
        value = rdata.conversion.convert(obj)
        if obj.info.type == RObjectType.STR and obj.attributes is not None:
            attributes = rdata.conversion.convert(obj.attributes)
            if "dim" in attributes:
                value = np.asarray(value, dtype=object).reshape(
                    tuple(attributes["dim"]), order="F"
                )
        return value

    return convert(parsed.object.value[0])


def format_matrix_column(values, digits=7):
    """Common precision used when R converts a numeric data-frame column to text."""
    from decimal import Decimal

    finite = [
        float(value)
        for value in values
        if value is not None and math.isfinite(float(value))
    ]
    if not finite:
        return ["NA"] * len(values)
    rounded = [format(value, f".{digits}g") for value in finite]
    fixed = [format(Decimal(value), "f") for value in rounded]
    decimals = max((len(value.split(".")[1]) if "." in value else 0) for value in fixed)
    significant = max(
        len(value.lstrip("-0").replace(".", "").lstrip("0").rstrip("0"))
        for value in fixed
    )
    significant = max(1, significant)
    as_fixed = [format(value, f".{decimals}f") for value in finite]
    as_scientific = [format(value, f".{significant - 1}e") for value in finite]
    scientific = max(map(len, as_scientific)) < max(map(len, as_fixed))
    output = []
    for value in values:
        if value is None or math.isnan(float(value)):
            output.append("NA")
        else:
            output.append(
                format(
                    float(value),
                    f".{significant - 1}e" if scientific else f".{decimals}f",
                )
            )
    width = max(map(len, output), default=0)
    return [value.rjust(width) for value in output]
