"""GEV energy fitting and IntaRNA p-values, translated from add_pval_to_csv_evdfit.R."""

from pathlib import Path
import math
import numpy as np

from .numerics import nelder_mead, r_mean, r_sum
from .statistics import _r_normalize_numeric_columns, r_number


def fit_gev(data):
    data = np.asarray(data, dtype=float)
    mean = r_mean(data)
    variance = r_sum((data - mean) ** 2) / (len(data) - 1)
    sigma = math.sqrt(6 * variance) / math.pi
    initial = [0.1, sigma, mean - 0.57722 * sigma]

    def objective(theta):
        xi, scale, location = theta
        with np.errstate(all="ignore"):
            y = 1 + xi * (data - location) / scale
            if scale < 0 or np.min(y) < 0:
                return 1e6
            return (
                len(data) * math.log(scale)
                + r_sum((1 + 1 / xi) * np.log(y))
                + r_sum(y ** (-1 / xi))
            )

    return nelder_mead(objective, initial)


def gev_pvalues(data, parameters):
    xi, scale, location = parameters
    with np.errstate(all="ignore"):
        return np.round(
            1
            - np.exp(-((1 + xi * (np.asarray(data) - location) / scale) ** (-1 / xi))),
            7,
        )


def add_energy_pvalues(directory: Path) -> str:
    options = (directory / "CopraRNA_option_file.txt").read_text().splitlines()
    region = options[3].split(":")[1]
    optimal = sorted(directory.glob("*_opt.intarna.csv"))
    suboptimal = sorted(directory.glob("*_subopt.intarna.csv"))
    for opt, sub in zip(optimal, suboptimal):
        if opt.name[:20] != sub.name[:20]:
            with (directory / "err.log").open("a") as stream:
                stream.write("Files don't match in add_pval_to_csv_evdfit.R\n")
        tables = []
        for path in (opt, sub):
            lines = path.read_text().splitlines()
            header = lines[0].split(";")
            rows = _r_normalize_numeric_columns([line.split(";") for line in lines[1:]])
            index = header.index("E_norm" if region == "cds" else "E")
            energies = [-float(row[index]) for row in rows]
            tables.append((path, header, rows, energies))
        parameters = fit_gev(tables[0][3])
        for path, header, rows, energies in tables:
            pvalues = gev_pvalues(energies, parameters)
            order = np.argsort(pvalues, kind="stable")
            header.append("p-value")
            path.write_text(
                ";".join(header)
                + "\n"
                + "".join(
                    ";".join(rows[i]) + ";" + r_number(float(pvalues[i])) + "\n"
                    for i in order
                )
            )
    return ""
