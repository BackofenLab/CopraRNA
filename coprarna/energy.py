"""GEV energy fitting and IntaRNA p-values, translated from add_pval_to_csv_evdfit.R."""

from pathlib import Path
from .io import paths, read_text, write_text
import math
import numpy as np

from .numerics import nelder_mead, r_mean, r_sum
from .statistics import _r_normalize_numeric_columns, r_number


def fit_gev(data):
    """Estimate (shape, scale, location) from negated interaction energies.

    Use the original moment-based initialization, likelihood penalty and R-style
    Nelder-Mead search. Replacing this with a generic SciPy fit changes captured
    scientific results; retained oracle fixtures constrain the implementation."""
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
    """Evaluate the fitted GEV upper tail and round to seven decimal places."""
    xi, scale, location = parameters
    with np.errstate(all="ignore"):
        return np.round(
            1
            - np.exp(-((1 + xi * (np.asarray(data) - location) / scale) ** (-1 / xi))),
            7,
        )


def add_energy_pvalues(directory: Path) -> str:
    """Fit optimal non-artifact energies and score both interaction sets.

    CDS predictions use length-normalized energy; other regions use raw energy.
    Read and rewrite gzip directly, retain the native schema plus p-value, and
    fail clearly when filtering leaves too little variable data to fit a GEV."""
    from .runtime import options

    region = options(directory)["region"]
    optimal = paths(directory, "*_opt.intarna.csv")
    suboptimal = paths(directory, "*_subopt.intarna.csv")
    for opt, sub in zip(optimal, suboptimal):
        if opt.name[:20] != sub.name[:20]:
            with (directory / "err.log").open("a") as stream:
                stream.write("Files don't match in add_pval_to_csv_evdfit.R\n")
        tables = []
        for path in (opt, sub):
            lines = read_text(path).splitlines()
            header = lines[0].split(";")
            rows = _r_normalize_numeric_columns([line.split(";") for line in lines[1:]])
            index = header.index("E_norm" if region == "cds" else "E")
            energies = [-float(row[index]) for row in rows]
            tables.append((path, header, rows, energies))
        if len(tables[0][3]) < 3 or len(set(tables[0][3])) < 2:
            raise ValueError(
                f"Insufficient distinct non-artifact interactions for GEV fitting in {opt.name}"
            )
        parameters = fit_gev(tables[0][3])
        for path, header, rows, energies in tables:
            pvalues = gev_pvalues(energies, parameters)
            order = np.argsort(pvalues, kind="stable")
            header.append("p-value")
            write_text(
                path,
                ";".join(header)
                + "\n"
                + "".join(
                    ";".join(rows[i]) + ";" + r_number(float(pvalues[i])) + "\n"
                    for i in order
                ),
            )
    return ""
