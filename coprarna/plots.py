"""Interaction-region and conservation plots using the legacy geometry and colors."""

import math
import os
from pathlib import Path
import re
import tempfile
import numpy as np
from Bio import SeqIO
from .runtime import options, command
from .rcompat import read_table, load_rdata
from .conservation import organism_names


def _pyplot():
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "coprarna-matplotlib")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


COLORS = [
    "#5cacee",
    "#ff4500",
    "#9acd32",
    "#b8860b",
    "#9a32cd",
    "#1874cd",
    "#cd661d",
    "#66cd00",
]


def _fasta(path):
    return {record.id: str(record.seq).lower() for record in SeqIO.parse(path, "fasta")}


def _positions(sequence):
    return np.flatnonzero(np.array(list(sequence)) != "-") + 1


def region_plot_data(path, number, directory):
    opt = options(directory)
    up = int(opt["nt upstream"])
    down = int(opt["nt downstream"])
    header, rows = read_table(path)
    end = next(i for i, name in enumerate(header) if "Annotation" in name)
    organisms = header[2:end]
    command(
        [
            "mafft",
            "--maxiterate",
            "1000",
            "--localpair",
            "--quiet",
            "--inputorder",
            "ncrna.fa",
        ],
        directory,
        output="ncrna_aligned.fa",
    )
    srna = {
        name.removeprefix("ncRNA_"): sequence
        for name, sequence in _fasta(directory / "ncrna_aligned.fa").items()
    }
    sequences = {
        org: _fasta(directory / f"{org}_upfromstartpos_{up}_down_{down}.fa")
        for org in organisms
    }
    result = {"mRNA": [], "sRNA": []}
    labels = []
    for i, row in enumerate(rows[:number]):
        entries = []
        name = ""
        for j, cell in enumerate(row[2:end]):
            tag = cell.split("(")[0] if cell else f"NA_{j + 1}"
            seq = next(
                (
                    value
                    for key, value in sequences[organisms[j]].items()
                    if re.search("i-" + tag, "i-" + key, re.I)
                ),
                "",
            )
            entries.append((tag, seq))
            if cell and not name:
                name = tag + "_" + cell.split("(", 1)[1].split("|")[0]
        labels.append(name)
        (directory / "fasta_temp_file").write_text(
            "".join(f">{tag}\n{seq}\n" for tag, seq in entries)
        )
        command(
            [
                "mafft",
                "--maxiterate",
                "1000",
                "--localpair",
                "--quiet",
                "--inputorder",
                "fasta_temp_file",
            ],
            directory,
            output="fasta_temp_file_out",
        )
        aligned = _fasta(directory / "fasta_temp_file_out")
        for j, (tag, _) in enumerate(entries):
            cell = row[j + 2]
            fields = cell.split("|") if cell else [""] * 8
            for feature, seq, startcol, endcol in [
                ("mRNA", aligned.get(tag, ""), 3, 4),
                ("sRNA", srna[organisms[j]], 5, 6),
            ]:
                positions = _positions(seq) if cell else np.array([0])
                gaps = (
                    np.flatnonzero(np.array(list(seq)) == "-") + 1
                    if cell
                    else np.array([0])
                )
                if not len(gaps):
                    gaps = np.array([0])
                start = int(float(fields[startcol])) if cell else 0
                stop = int(float(fields[endcol])) if cell else 0
                interaction = (
                    [int(positions[start - 1]), int(positions[stop - 1])]
                    if cell
                    else [0, 0]
                )
                shift = (
                    up * 1.2 - int(positions[up - 1])
                    if feature == "mRNA" and cell
                    else 0
                )
                result[feature].append(
                    dict(
                        start=float(min(positions) + shift),
                        end=float(max(positions) + shift),
                        interaction=[value + shift for value in interaction],
                        gaps=(gaps + shift).tolist(),
                        y=j + 1 + i * (len(organisms) + 1),
                        cluster=i,
                    )
                )
    return result, labels, organisms


def plot_interaction_regions(path: Path, number: int, directory: Path) -> str:
    opt = options(directory)
    # Retain the legacy lowercase-CDS plotting failure; fixed separately.
    from .compatibility import legacy_enabled

    if opt["region"] == "cds" and legacy_enabled():
        raise ValueError(
            "Legacy region plotting does not initialize center for lowercase cds"
        )
    header, rows = read_table(path)
    if len(rows) < max(15, number):
        raise ValueError("Legacy region plotting indexes beyond available predictions")
    number = max(15, number)
    data, labels, orgs = region_plot_data(path, number, directory)
    plt = _pyplot()
    up = int(opt["nt upstream"])
    size = len(orgs)
    from matplotlib.backends.backend_pdf import PdfPages

    def draw(feature, points, first, count, histogram):
        fig, axis = plt.subplots(figsize=(11.69, 8.27))
        mi = min(p["start"] for p in points)
        ma = max(p["end"] for p in points)
        for point in points:
            y = point["y"] - first * (size + 1)
            axis.plot(
                [point["start"], point["end"]], [y, y], color="#d9d9d9", linewidth=1.5
            )
            axis.plot(
                point["interaction"],
                [y, y],
                color=COLORS[point["cluster"] % 8],
                linewidth=1.5,
            )
            for gap in point["gaps"]:
                axis.plot([gap, gap + 1], [y, y], color="white", linewidth=1.5)
        for i in range(count):
            y = (i + 1) * (size + 1)
            axis.plot([mi, ma], [y, y], color="black", linewidth=0.5)
            axis.text(
                ma * 1.15, y - size / 2, labels[first + i], fontsize=8, ha="center"
            )
        limit = count * size + count - 1
        axis.set_xlim(mi, ma * 1.25)
        axis.set_ylim(0, limit * (2.02 if histogram else 1.01))
        axis.set_yticks([])
        if histogram:
            positions = []
            for point in data[feature][: min(number, 100) * size]:
                start, end = point["interaction"]
                gaps = set(point["gaps"])
                positions.extend(
                    value for value in np.arange(start, end + 1) if value not in gaps
                )
            unique, counts = np.unique(positions, return_counts=True)
            if len(unique):
                xs = np.arange(min(unique), max(unique) + 1)
                ys = np.zeros(len(xs))
                ys[(unique - min(unique)).astype(int)] = counts
                base = limit * 1.05
                ys = ys / max(ys) * limit + base
                axis.fill_between(
                    np.r_[xs[0] - 1, xs, xs[-1] + 1],
                    base,
                    np.r_[base, ys, base],
                    facecolor="white",
                    edgecolor="black",
                    linewidth=0.7,
                )
                axis.axhline(base, color="black", linewidth=0.5)
        if feature == "mRNA" and opt["region"] != "cds":
            center = up * 1.2
            axis.axvline(center, color="gray", ymax=0.49 if histogram else 1)
            offsets = (
                [-200, -100, -50, 0, 50, 100, 200]
                if opt["region"] == "5utr"
                else [-100, 0, 100]
            )
            ticklabels = (
                [str(x) if x else "start codon" for x in offsets]
                if opt["region"] == "5utr"
                else ["100", "stop codon", "100"]
            )
            axis.set_xticks([center + x for x in offsets], ticklabels)
        else:
            axis.set_xlabel("nt")
        axis.set_title(
            "Predicted interaction regions "
            + ("mRNAs" if feature == "mRNA" else "sRNA")
        )
        for side in ("top", "left", "right"):
            axis.spines[side].set_visible(False)
        fig.tight_layout()
        return fig

    for feature in ("mRNA", "sRNA"):
        prefix = feature + "_regions_with_histogram"
        fig = draw(feature, data[feature][: 15 * size], 0, 15, True)
        for extension in ("pdf", "ps", "png"):
            fig.savefig(directory / f"{prefix}.{extension}", dpi=300)
        if int(opt["websrv"]):
            fig.savefig(directory / f"thumbnail_{feature}.png", dpi=20)
        plt.close(fig)
        prefix = (
            "mRNA_regions_extended_list"
            if feature == "mRNA"
            else "sRNA_regions_with_extended_list"
        )
        with PdfPages(directory / f"{prefix}.pdf") as pdf:
            for first in range(0, number, 25):
                count = min(25, number - first)
                fig = draw(
                    feature,
                    data[feature][first * size : (first + count) * size],
                    first,
                    count,
                    False,
                )
                pdf.savefig(fig)
                plt.close(fig)
        command(
            [
                "gs",
                "-q",
                "-dBATCH",
                "-dNOPAUSE",
                "-sDEVICE=ps2write",
                f"-sOutputFile={prefix}.ps",
                f"{prefix}.pdf",
            ],
            directory,
        )
    return ""


def evolutionary_heatmap_data(directory: Path, **settings):
    opt = options(directory)
    inputfile = settings.get("inputfile", "CopraRNA_result_all.csv")
    number = int(settings.get("num", 25))
    header, rows = read_table(directory / inputfile)
    end = next(i for i, name in enumerate(header) if "Annotation" in name)
    orgs = header[2:end]
    selected = rows[:number]
    if str(settings.get("select", "FALSE")).upper() == "TRUE":
        genes = (
            (directory / settings.get("sel", "genelist.txt")).read_text().splitlines()
        )
        selected = [
            row
            for row in rows
            if any(re.search(gene, "".join(row[2:end])) for gene in genes)
        ]
    tables = load_rdata(directory / "conservation_table.Rdata", "conservation_table")
    ooi = int(opt["cons"]) == 1 or settings.get("consensus") == "ooi"
    optimal, suboptimal = tables[2:4] if ooi else tables[:2]
    scores = np.full((len(selected), len(orgs)), np.nan)
    notes = np.full(scores.shape, "", dtype=object)
    for i, row in enumerate(selected):
        index = int(float(row[-1])) - 1
        for j in range(len(orgs)):

            def decode(value):
                if value is None:
                    return None, math.nan
                fields = str(value).split("|")
                if len(fields) != 4:
                    return None, math.nan
                return (1 if fields[1] == "TRUE" else 0) + (
                    2 if "TRUE" in fields[2:] else 0
                ), float(fields[0])

            code, p = decode(optimal[index, j])
            subcode, sp = decode(suboptimal[index, j])
            if code in (0, 2) and subcode in (1, 3):
                code = 4 if subcode == 1 else 5
                p = sp
            if code is not None:
                scores[i, j] = (
                    np.searchsorted(
                        [
                            0.001,
                            0.01,
                            0.1,
                            0.2,
                            float(settings.get("int_p_thres", 0.3)),
                        ],
                        p,
                        side="left",
                    )
                    + 1
                )
                notes[i, j] = {0: "-/-", 1: "+/-", 2: "-/+", 3: "", 4: "2+/-", 5: "2"}[
                    code
                ]
    labels = organism_names(orgs, directory)
    clustering = settings.get("clustering", "sRNA")
    tree = None
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform, pdist

    if clustering == "sRNA":
        command(
            [
                "clustalo",
                "-i",
                "input_sRNA.fa",
                "--distmat-out=distmatout1.txt",
                "--full",
                "--use-kimura",
                "--output-order=input-order",
                "--force",
                "--max-hmm-iterations=-1",
            ],
            directory,
        )
        lines = (directory / "distmatout1.txt").read_text().splitlines()[1:]
        identifiers = [line.split()[0] for line in lines]
        distances = np.array([[float(v) for v in line.split()[1:]] for line in lines])
        (directory / "distmatout1.txt").unlink()
        indices = [
            next(i for i, org in enumerate(orgs) if org.split(".")[0] == name)
            for name in identifiers
        ]
        scores = scores[:, indices]
        notes = notes[:, indices]
        tree = linkage(squareform(distances, checks=False), method="average")
    elif clustering == "ribosomal":
        lines = (directory / "distmat.out").read_text().splitlines()[7:]
        cells = [line.split("\t") for line in lines if line.strip()]
        identifiers = [line[-1].split()[0] for line in cells]
        matrix = np.array(
            [[float(v) if v.strip() else 0 for v in row[1:-2]] for row in cells]
        )
        distances = matrix + matrix.T
        indices = [
            next(i for i, org in enumerate(orgs) if org.split(".")[0] == name)
            for name in identifiers
        ]
        scores = scores[:, indices]
        notes = notes[:, indices]
        tree = linkage(squareform(distances, checks=False), method="average")
    else:
        tree = linkage(pdist(np.nan_to_num(scores.T)), method="complete")
    genes = []
    for row in selected:
        cell = next((cell for cell in row[2:end] if cell), "")
        genes.append(
            cell.split("(", 1)[-1].split("|")[0].replace("/", "")
            + "_"
            + cell.split("(")[0]
        )
    # Match R hclust: put the singleton/tighter subtree on the left.
    for i in range(len(tree)):
        a, b = map(int, tree[i, :2])
        n = len(orgs)

        def key(node):
            return (-1, node) if node < n else (tree[node - n, 2], node)

        if key(a) > key(b):
            tree[i, 0], tree[i, 1] = b, a
    return scores, notes, tree, [a + "_" + b for a, b in zip(labels, orgs)], genes


def plot_evolutionary_heatmap(directory: Path, **settings) -> str:
    scores, notes, tree, labels, genes = evolutionary_heatmap_data(
        directory, **settings
    )
    plt = _pyplot()
    from matplotlib.colors import ListedColormap
    from scipy.cluster.hierarchy import dendrogram, leaves_list

    order = leaves_list(tree)
    scores = scores[:, order]
    notes = notes[:, order]
    palette = ListedColormap(
        ["#caff70", "#c0ff3e", "#b3ee3a", "#9acd32", "#698b22", "#ff8c69"]
    )
    palette.set_bad("black")
    fig = plt.figure(figsize=(11.69, 8.27))
    layout = fig.add_gridspec(2, 1, height_ratios=[1, 8])
    top = fig.add_subplot(layout[0])
    dendrogram(
        tree, ax=top, no_labels=True, color_threshold=0, above_threshold_color="black"
    )
    top.axis("off")
    axis = fig.add_subplot(layout[1])
    image = axis.imshow(scores, cmap=palette, vmin=0.5, vmax=6.5, aspect="auto")
    for i in range(len(genes)):
        for j in range(len(labels)):
            axis.text(
                j, i, notes[i, j], color="red", ha="center", va="center", fontsize=6
            )
    axis.set_xticks(
        range(len(labels)),
        [labels[i] for i in order],
        rotation=45,
        ha="right",
        fontsize=8,
    )
    axis.set_yticks(range(len(genes)), genes, fontsize=7)
    axis.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.5)
    axis.tick_params(which="minor", length=0)
    mode = (
        "ooi"
        if int(options(directory)["cons"]) == 1
        else settings.get("consensus", "overall")
    )
    top.set_title(mode + " _consensus_based")
    legend = fig.colorbar(image, ax=axis, ticks=range(1, 7), fraction=0.025, pad=0.02)
    threshold = settings.get("int_p_thres", 0.3)
    legend.ax.set_yticklabels(
        ["<0.001", "<0.01", "<0.1", "<0.2", f"<{threshold}", f">{threshold}"]
    )
    fig.tight_layout()
    fig.savefig(
        directory / (settings.get("prefix", "sRNA") + "_conservation_heatmap.pdf")
    )
    plt.close(fig)
    return ""
