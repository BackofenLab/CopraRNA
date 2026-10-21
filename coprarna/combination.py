"""CopraRNA p-value combination and organism subset selection."""

from __future__ import annotations
from collections import OrderedDict
from pathlib import Path
import math
import numpy as np
from scipy import special
from scipy.linalg.blas import dnrm2, ddot, daxpy

from .numerics import r_sum, normal_upper, normal_probit
from .rcompat import read_table, write_table, numeric_columns, make_names, load_rdata
from .statistics import bh_adjust


def fit_rho(probits, weights):
    """Fit inter-organism correlation against a uniform empirical p-value CDF.

    The one-parameter least-squares iteration follows R nls: finite-difference
    gradient, LINPACK-style Householder rounding and step halving. BLAS reduction
    order matters for byte parity; do not replace this with a generic optimizer."""
    a = np.array([r_sum(row * weights) for row in probits])
    d = r_sum(np.asarray(weights) ** 2)
    b = r_sum(weights) ** 2
    y = (1 / len(a)) * np.arange(1, len(a) + 1, dtype=float)

    def prediction(rho):
        variance = (1 - rho) * d + rho * b
        if variance <= 0:
            return np.full(len(a), np.nan)
        return np.sort(normal_upper(a / math.sqrt(variance)))

    rho = 0.0
    factor = 1.0
    predicted = prediction(rho)
    residual = y - predicted
    deviance = r_sum(residual**2)
    for _ in range(50):
        step = math.sqrt(np.finfo(float).eps) * (abs(rho) or 1.0)
        gradient = (prediction(rho + step) - predicted) / step
        norm = r_sum(gradient**2)
        if not np.isfinite(norm) or norm == 0:
            break
        # nls uses a one-column LINPACK QR decomposition, including its
        # Householder rounding, rather than a normal-equation dot product.
        magnitude = float(dnrm2(gradient))
        signed = math.copysign(magnitude, float(gradient[0]))
        householder = gradient * (1 / signed)
        householder[0] += 1
        transformed = daxpy(
            householder,
            residual.copy(),
            a=-float(ddot(householder, residual)) / householder[0],
        )
        delta = -float(transformed[0]) / signed
        denominator = r_sum(transformed[1:] ** 2)
        if denominator > 0 and abs(transformed[0]) / math.sqrt(denominator) <= 1e-5:
            break
        accepted = False
        while factor >= 1 / 128:
            proposed = rho + factor * delta
            fitted = prediction(proposed)
            next_deviance = r_sum((y - fitted) ** 2)
            if np.isfinite(next_deviance) and next_deviance <= deviance:
                rho, predicted, deviance = proposed, fitted, next_deviance
                residual = y - predicted
                factor = min(1.0, 2 * factor)
                accepted = True
                break
            factor /= 2
        if not accepted:
            break
    return rho


def _combine_matrix(matrix, weights, selection=None):
    """Combine p-values for all eligible organism subsets of size at least three.

    Matrix rows are target clusters and columns are organisms. Selection may use
    conserved interactions; fit correlation on complete rows for each subset and
    keep the minimum combined probability for each cluster and organism."""
    matrix = np.asarray(matrix, float)
    selected = matrix if selection is None else np.asarray(selection, float)
    variants = OrderedDict()
    for i, row in enumerate(selected):
        indices = np.argsort(row, kind="stable")
        indices = indices[np.isfinite(row[indices])]
        for size in range(3, len(indices) + 1):
            key = tuple(sorted(indices[:size]))
            variants.setdefault(key, []).append(i)
    out = np.full(matrix.shape, np.nan)
    balanced = np.full(len(matrix), np.nan)
    for indices, positions in variants.items():
        indices = list(indices)
        full = matrix[:, indices]
        valid = np.flatnonzero(~np.isnan(full).any(axis=1))
        probits = normal_probit(full[valid])
        w = weights[indices]
        rho = fit_rho(probits, w)
        d = r_sum(w * w)
        b = r_sum(w) ** 2
        denominator = (1 - rho) * d + rho * b
        source = selected if selection is not None else matrix
        for position in positions:
            p = source[position, indices]
            with np.errstate(all="ignore"):
                value = float(
                    normal_upper(r_sum(w * normal_probit(p)) / math.sqrt(denominator))
                )
            balanced[position] = np.fmin(balanced[position], value)
            out[position, indices] = np.fmin(out[position, indices], value)
    return balanced, out


def _rra(rank_lists):
    """Aggregate ranked cluster lists using corrected beta order-statistic scores.

    Missing entries receive rank one. The small-argument beta expansion preserves
    the numerical convention of the captured RobustRankAggreg output."""
    universe = list(
        dict.fromkeys(index for sequence in rank_lists for index in sequence)
    )
    total = len(universe)
    n = len(rank_lists)
    scores = {}
    for index in universe:
        ranks = sorted(
            (sequence.index(index) + 1) / total if index in sequence else 1.0
            for sequence in rank_lists
        )

        def beta(i, value):
            a = float(i)
            b = float(n - i + 1)
            if value == 0:
                return 0.0
            if value == 1:
                return 1.0
            if b * value <= 0.7 or b <= 1:
                front = math.exp(a * math.log(value) - float(special.betaln(a, b))) / a
                coefficient = 1.0
                total = 0.0
                for k in range(1, n - i + 2):
                    coefficient *= (0.5 - b / k + 0.5) * value
                    total += coefficient / (a + k)
                    if coefficient == 0:
                        break
                return front * (1 + a * total)
            return float(special.betainc(i, n - i + 1, value))

        scores[index] = min(
            1.0, n * min(beta(i + 1, value) for i, value in enumerate(ranks))
        )
    return sorted(universe, key=lambda index: scores[index]), scores


def _conservation_p(table, subtable, require_srna):
    """Choose the first optimal/suboptimal probability satisfying consensus flags."""
    result = np.full(np.shape(table), np.nan)
    for i in range(len(table)):
        for j in range(len(table[i])):
            for data in (table, subtable):
                value = data[i][j]
                if value is None:
                    continue
                fields = str(value).split("|")
                if (
                    len(fields) == 4
                    and fields[1] == "TRUE"
                    and (not require_srna or "TRUE" in fields[2:])
                ):
                    result[i, j] = float(fields[0])
                    break
    return result


def join_copra2(ooi: str, directory: Path, *modes: str) -> str:
    """Write balanced, organism-of-interest and consensus comparative rankings.

    Read annotated interactions, phylogenetic weights and persisted conservation
    matrices. Apply the configured root to weights, combine eligible subsets,
    adjust p-values and preserve the report column/ordering conventions."""
    header, raw = read_table(
        directory / "CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv"
    )
    header = make_names(header)
    rows = numeric_columns(raw)
    end = next(i for i, name in enumerate(header) if "Annotation" in name)
    organisms = header[2:end]
    matrix = np.full((len(rows), len(organisms)), np.nan)
    for i, row in enumerate(raw):
        for j, cell in enumerate(row[2:end]):
            if cell:
                matrix[i, j] = float(cell.split("|")[2])
    from .runtime import options

    root = float(options(directory)["root"])
    _, weight_rows = read_table(directory / "zscore.weight", sep=";")
    # Weight files have no header.
    weight_rows = [
        line.split(";")
        for line in (directory / "zscore.weight").read_text().splitlines()
    ]
    mapping = {name.upper(): float(w) for name, w in weight_rows}
    weights = np.array([mapping[org] for org in organisms])
    weights = np.ones(len(weights)) if root == 0 else weights ** (1 / root)
    balanced, per_organism = _combine_matrix(matrix, weights)
    focus = next(i for i, name in enumerate(organisms) if ooi in name)
    output_header = ["fdr", "p-value"] + header[2:] + ["initial_sorting"]

    def write(name, values):
        adjusted = bh_adjust(list(values))
        order = sorted(
            range(len(values)),
            key=lambda i: (
                math.isnan(values[i]),
                values[i] if not math.isnan(values[i]) else 0,
            ),
        )
        write_table(
            directory / name,
            output_header,
            [[adjusted[i], values[i], *rows[i][2:], i + 1] for i in order],
        )

    write("CopraRNA2_final_all_balanced.csv", balanced)
    write("CopraRNA2_final_all_ooi.csv", per_organism[:, focus])
    ranks = np.full(matrix.shape, np.nan)
    top = []
    for j in range(len(organisms)):
        valid = np.flatnonzero(~np.isnan(per_organism[:, j]))
        order = valid[np.argsort(per_organism[valid, j], kind="stable")]
        ranks[order, j] = np.arange(1, len(order) + 1)
        top.extend(int(i) for i in order[:200])
    top = list(dict.fromkeys(top))
    rank_lists = [
        sorted(
            [k for k, i in enumerate(top) if np.isfinite(ranks[i, j])],
            key=lambda k: ranks[top[k], j],
        )
        for j in range(len(organisms))
    ]
    rank_order, rank_scores = _rra(rank_lists)
    evo_header = (
        ["Rank_p_value", "evo_ooi_rank_result", "evo_ooi_p_result", "evo_int_p_result"]
        + header[2:]
        + organisms * 3
        + ["initial_sorting"]
    )
    evo_rows = []
    for k in rank_order:
        i = top[k]
        evo_rows.append(
            [
                rank_scores[k],
                float(np.sum(ranks[i] <= 250)) / len(organisms),
                float(np.sum(per_organism[i] <= 0.001)) / len(organisms),
                float(np.sum(matrix[i] <= 0.3)) / len(organisms),
                *rows[i][2:],
                *ranks[i],
                *per_organism[i],
                *matrix[i],
                i + 1,
            ]
        )
    write_table(
        directory / "CopraRNA2_final_all_evo.csv",
        evo_header,
        evo_rows,
        sep="\t",
        quote=True,
    )
    if "overall_consensus" in modes or "ooi_consensus" in modes:
        tables = load_rdata(
            directory / "conservation_table.Rdata", "conservation_table"
        )
        tables = [np.asarray(table, dtype=object) for table in tables]
        if "overall_consensus" in modes:
            selection = _conservation_p(tables[0], tables[1], "sRNA" in modes)
            balanced, per = _combine_matrix(matrix, weights, selection)
            write("CopraRNA2_final_all_balanced_consensus.csv", balanced)
            write("CopraRNA2_final_all_ooi_consensus.csv", per[:, focus])
        if "ooi_consensus" in modes:
            selection = _conservation_p(tables[2], tables[3], "sRNA" in modes)
            _, per = _combine_matrix(matrix, weights, selection)
            write("CopraRNA2_final_all_ooi_ooiconsensus.csv", per[:, focus])
    return ""


def combine_clusters(organism_count: int, directory: Path) -> str:
    """Build the 16S tree, SVG and weights, then prepare comparative cluster rows.

    Temporarily shorten long accessions for PHYLIP's name limit and restore them
    in both tree/report outputs before calculating organism weights."""
    from .runtime import command
    from .phylogeny import compute_weights, render_tree
    from .statistics import prepare_copra2
    from .intarna import numeric_sort

    text = (directory / "compatible.distmat").read_text()
    mapped = text
    mapping = {}
    for line in text.splitlines()[1:]:
        identifier = line.split("\t")[0]
        if len(identifier) > 10:
            short = f"na_{1000 + len(mapping)}"
            mapping[short] = identifier
            mapped = mapped.replace(identifier, short)
    (directory / "compatible.distmat.mapped").write_text(mapped)
    command(
        [
            "fneighbor",
            "-datafile",
            "compatible.distmat.mapped",
            "-outfile",
            "compatible.fneighbor.mapped",
            "-treetype",
            "u",
            "-auto",
        ],
        directory,
    )
    for source, dest in [
        ("compatible.fneighbor.mapped", "compatible.fneighbor"),
        ("distmat.treefile", "compatible.treefile"),
    ]:
        text = (directory / source).read_text().replace("0.00000", "0.00001")
        for short, identifier in mapping.items():
            text = text.replace(short, identifier)
        (directory / dest).write_text(text)
        (directory / source).unlink()
    render_tree(
        directory / "compatible.treefile", directory / "compatible.treefile.svg"
    )
    (directory / "zscore.weight").write_text(
        compute_weights(
            directory / "compatible.treefile", directory / "compatible.fneighbor"
        )
    )
    prepare_copra2(directory / "opt_tags.clustered", directory)
    path = directory / "CopraRNA2_prep.csv"
    (directory / "CopraRNA2_prep_sorted.csv").write_text(
        "".join(numeric_sort(path.read_text().splitlines(keepends=True)))
    )
    return ""
