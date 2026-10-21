"""R table serialization and persisted matrix conventions used by CopraRNA."""

from pathlib import Path
import csv
import math
import re
import numpy as np
from .statistics import r_number


def read_table(path, sep=","):
    with Path(path).open(newline="") as stream:
        rows = list(csv.reader(stream, delimiter=sep))
    return rows[0], rows[1:]


def make_names(names):
    return [re.sub(r"[^\w.]", ".", name) for name in names]


def scalar(value, quote=False):
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
    lines = []
    if header is not None:
        lines.append(
            sep.join(
                scalar(v, quote) for v in ([""] if row_names else []) + list(header)
            )
        )
    for i, row in enumerate(rows, 1):
        cells = ([str(i)] if row_names else []) + list(row)
        lines.append(sep.join(scalar(v, quote) for v in cells))
    Path(path).write_text("\n".join(lines) + "\n")


def numeric_columns(rows):
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
    def __init__(self, values, columns=None, rows=None):
        self.values = np.asarray(values, dtype=object)
        self.columns = columns
        self.rows = rows


def save_rdata(path, name, values, *, names=None, columns=None):
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

    if columns is not None:
        values = [RMatrix(value, columns) for value in values]
    rdata.write_rda(path, {name: values}, constructor_dict={RMatrix: constructor})


def load_rdata(path, name):
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
