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
from .statistics import bh_adjust, r_number


def fit_rho(probits, weights):
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
    root = float(
        (directory / "CopraRNA_option_file.txt")
        .read_text()
        .splitlines()[13]
        .split(":")[1]
    )
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


class RRandom:
    """Mersenne Twister initialization and inversion-normal sampling used by R."""

    def __init__(self, seed=None):
        if seed is None:
            import secrets

            seed = secrets.randbits(32)
        seed = int(seed) & 0xFFFFFFFF
        for _ in range(50):
            seed = (69069 * seed + 1) & 0xFFFFFFFF
        state = []
        for _ in range(625):
            seed = (69069 * seed + 1) & 0xFFFFFFFF
            state.append(seed)
        self.random = np.random.RandomState()
        self.random.set_state(
            ("MT19937", np.array(state[1:], dtype=np.uint32), 624, 0, 0.0)
        )

    def uniform(self):
        value = int(self.random.randint(0, 2**32, dtype=np.uint32)) / 2**32
        return max(0.5 / 2**32, min(1 - 0.5 / 2**32, value))

    def normal(self):
        value = (math.floor(134217728 * self.uniform()) + self.uniform()) / 134217728
        return -float(normal_probit(value))


def _histogram_error(values):
    values = np.asarray(values, float)
    lo, hi = float(np.min(values)), float(np.max(values))
    cell = (hi - lo) / 100
    if not cell:
        return math.inf
    base = 10 ** math.floor(math.log10(cell))
    unit = base
    if 2 * base - cell < 1.5 * (cell - base):
        unit = 2 * base
    if 5 * base - cell < 2.75 * (cell - unit):
        unit = 5 * base
    if 10 * base - cell < 1.5 * (cell - unit):
        unit = 10 * base
    breaks = np.arange(math.floor(lo / unit), math.ceil(hi / unit) + 1) * unit
    fuzz = 1e-7 * float(np.median(np.diff(breaks)))
    fuzzy = breaks + fuzz
    fuzzy[0] -= 2 * fuzz
    bins = np.searchsorted(fuzzy, values, side="left") - 1
    _counts = np.bincount(np.clip(bins, 0, len(breaks) - 2), minlength=len(breaks) - 1)
    # Legacy hist objects have density, not intensities. Accessing the absent
    # intensities member yields NULL, and sum((1-NULL)^2) is zero.
    from .compatibility import legacy_enabled

    if legacy_enabled():
        return 0.0
    density = _counts / (len(values) * np.diff(breaks))
    return r_sum((1 - density) ** 2)


def join_copra1(path: Path, directory: Path, seed=None) -> str:
    _, rows = read_table(path, sep=";")
    options = (directory / "CopraRNA_option_file.txt").read_text().splitlines()
    root = float(options[13].split(":")[1])
    orgs = sorted(
        [
            line.split(";")[0]
            for line in (directory / "zscore.weight").read_text().splitlines()
        ]
    )
    wmap = dict(
        line.split(";")
        for line in (directory / "zscore.weight").read_text().splitlines()
    )
    weights = np.array([float(wmap[name]) ** (1 / root) for name in orgs])
    groups = OrderedDict()
    for row in rows:
        groups.setdefault(row[-1], []).append(row)
    full = []
    for group in groups.values():
        if len(group) == len(orgs):
            group = sorted(group, key=lambda row: row[1])
            full.append([normal_probit(float(row[35])) for row in group])
    full = np.asarray(full)
    full = full[np.isfinite(full.sum(axis=1))]
    means = full.mean(axis=0)
    covariance = np.cov(full, rowvar=False)
    random = RRandom(seed)

    def search(interval):
        best_error = 1e9
        best_rho = float(interval[0])
        best_lines = []
        for rho in interval:
            values = []
            lines = []
            for group in groups.values():
                group = sorted(group, key=lambda row: row[1])
                present = [orgs.index(row[1].replace("ncrna_", "")) for row in group]
                probabilities = np.array([float(row[35]) for row in group])
                probabilities[probabilities == 1] = 0.9999999999999999
                p = np.full(len(orgs), np.nan)
                p[present] = normal_probit(probabilities)
                missing = [j for j in range(len(orgs)) if j not in present]
                if missing:
                    c = covariance[np.ix_(missing, present)] @ np.linalg.inv(
                        covariance[np.ix_(present, present)]
                    )
                    conditional_mean = means[missing] + c @ (
                        p[present] - means[present]
                    )
                    conditional_cov = (
                        covariance[np.ix_(missing, missing)]
                        - c @ covariance[np.ix_(present, missing)]
                    )
                    for k, j in enumerate(missing):
                        p[j] = (
                            conditional_mean[k]
                            + math.sqrt(conditional_cov[k, k]) * random.normal()
                        )
                with np.errstate(all="ignore"):
                    value = float(
                        normal_upper(
                            r_sum(weights * p)
                            / np.sqrt(
                                (1 - rho) * r_sum(weights**2)
                                + rho * r_sum(weights) ** 2
                            )
                        )
                    )
                if math.isnan(value):
                    value = 0.0
                values.append(value)
                lines.append(
                    r_number(value) + ";" + ";".join(row[0] for row in group) + ";\n"
                )
            error = _histogram_error(values)
            if error < best_error:
                best_error, best_rho, best_lines = error, float(rho), lines
        return best_rho, best_lines

    rho, _ = search(np.arange(11) / 10)
    if rho == 0:
        rho = 0.1
    _, lines = search(np.arange(21) / 100 + rho - 0.1)
    (directory / "CopraRNA1_with_pvsample.csv").write_text("".join(lines))
    return ""


def combine_clusters(organism_count: int, directory: Path, seed=None) -> str:
    from .runtime import command, options
    from .phylogeny import compute_weights
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
    (directory / "zscore.weight").write_text(
        compute_weights(
            directory / "compatible.treefile", directory / "compatible.fneighbor"
        )
    )
    if int(options(directory)["CopraRNA1"]):
        join_copra1(directory / "opt_tags.clustered_rcsize", directory, seed)
        path = directory / "CopraRNA1_with_pvsample.csv"
        (directory / "CopraRNA1_with_pvsample_sorted.csv").write_text(
            "".join(numeric_sort(path.read_text().splitlines(keepends=True)))
        )
    prepare_copra2(directory / "opt_tags.clustered", directory)
    path = directory / "CopraRNA2_prep.csv"
    (directory / "CopraRNA2_prep_sorted.csv").write_text(
        "".join(numeric_sort(path.read_text().splitlines(keepends=True)))
    )
    return ""
