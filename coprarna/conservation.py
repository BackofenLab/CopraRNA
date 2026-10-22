"""Interaction-coordinate mapping, consensus peaks and Jalview annotations."""

from pathlib import Path
from .io import paths, genome_paths, genbank_records
import json
import numpy as np
from Bio import SeqIO

from .runtime import command, options, srna_input
from .rcompat import read_table, write_table, save_rdata, format_matrix_column
from .rcompat import RNamedList
from .statistics import r_number

COLORS = [
    "p<=0001\tD9F69C",
    "p<=001\tC2F161",
    "p<01\tBAEF4D",
    "p<02\tB3EE3A",
    "p<03\tA1D634",
    "p>03\tFF4500",
]
LABELS = ["p<=0001", "p<=001", "p<=01", "p<=02", "p<=03", "p>03"]
FIELDS = [
    "start",
    "end",
    "start_seed",
    "end_seed",
    "start_sRNA",
    "end_sRNA",
    "start_seed_sRNA",
    "end_seed_sRNA",
    "name",
    "name2",
    "pvalue",
    "orgs",
]
COORDS = [8, 9, 25, 26, 10, 11, 27, 28]


def peak_find(values, threshold=0.4):
    """Return inclusive 1-based boundaries around all maximum-height positions."""
    values = np.asarray(values)
    maximum = max(values)
    cutoff = maximum * threshold
    maxima = np.flatnonzero(values == maximum)
    start, end = int(maxima[0]), int(maxima[-1])
    while start >= 0 and values[start] >= cutoff:
        start -= 1
    while end < len(values) and values[end] >= cutoff:
        end += 1
    return [start + 2, end]


def peak_find_two(values, threshold=0.4, second_threshold=0.9):
    """Find the main consensus peak and an eligible near-height secondary peak.

    Retain the original asymmetric end threshold for reproducible peak boundaries."""
    first = peak_find(values, threshold)
    second = []
    outside = np.r_[values[: first[0] - 1], values[first[1] :]]
    unique = np.unique(outside)
    if len(unique) > 10:
        value = max(unique)
        maxima = np.flatnonzero(values == value)
        if (
            not any(first[0] - 1 <= i < first[1] for i in maxima)
            and value / max(values) >= second_threshold
        ):
            start, end = int(maxima[0]), int(maxima[-1])
            while start >= 0 and values[start] >= threshold * value:
                start -= 1
            # The original uses the first peak's threshold for the second end.
            while end < len(values) and values[end] >= threshold * max(values):
                end += 1
            second = [start + 2, end]
    return [first, second]


def overlap(peak, start, end, threshold=0.6):
    """Test overlap relative to the shorter inclusive interval; an absent peak fails."""
    if not peak:
        return False
    left, right = min(start, end), max(start, end)
    width = max(0, min(right, peak[1]) - max(left, peak[0]) + 1)
    return width / min(right - left + 1, abs(peak[1] - peak[0]) + 1) >= threshold


def organism_names(organisms, directory):
    """Return compact plot labels from the job metadata or local GenBank records.

    Match complete accession identities, ignoring version only for display lookup.
    Standalone reports fall back to the accession when metadata is unavailable.
    Read files without querying NCBI or modifying the job.
    """
    mapping = {}
    path = directory / "genome_metadata.json"
    if path.is_file():
        for item in json.loads(path.read_text()):
            for identifier in item["replicons"]:
                mapping[identifier.split(".")[0]] = item["organism"]
    else:
        for path in genome_paths(directory):
            for record in genbank_records(path):
                mapping[record.id.split(".")[0]] = record.annotations.get(
                    "organism", record.id
                )
    names = []
    for org in organisms:
        found = mapping.get(org.split(".")[0])
        if not found:
            names.append(org)
            continue
        parts = found.replace("_", " ").split()
        label = parts[0][:3] + "_" + (parts[1] if len(parts) > 1 else "NA")
        if len(parts) > 2:
            label += "_" + parts[-1]
        names.append(label)
    return names


def _fasta(path):
    return {record.id: str(record.seq).lower() for record in SeqIO.parse(path, "fasta")}


def _write_fasta(path, records):
    path.write_text(
        "".join(
            ">"
            + name
            + "\n"
            + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
            + "\n"
            for name, seq in records
        )
    )


def _boolean(value):
    return "TRUE" if value else "FALSE"


def compute_conservation(ooi: str, directory: Path) -> str:
    """Map interactions to target/sRNA alignments and persist consensus annotations.

    Combine optimal and suboptimal coordinates, locate weighted consensus peaks,
    and write R-compatible matrices plus Jalview annotation/alignment files. The
    OOI selects organism-specific mappings; other outputs retain organism order."""
    command(
        [
            "mafft",
            "--maxiterate",
            "1000",
            "--localpair",
            "--quiet",
            "--localpair",
            "--quiet",
            srna_input(directory),
        ],
        directory,
        output="aligned_sRNA.fa",
    )
    srna = _fasta(directory / "aligned_sRNA.fa")
    header, rows = read_table(
        directory / "CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv"
    )
    end = next(i for i, name in enumerate(header) if "Annotation" in name)
    organisms = header[2:end]
    labels = organism_names(organisms, directory)
    opt_header, opt_rows = read_table(directory / "opt_tags.clustered", sep=";")
    optimal_indices = {
        row[0].lower(): i for i, row in reversed(list(enumerate(opt_rows, 1)))
    }
    optimal = {}
    for row in opt_rows:
        optimal.setdefault(row[0].lower(), row)
    suboptimal = {}
    for path in paths(directory, "*_subopt.intarna.csv"):
        _, data = read_table(path, sep=";")
        for row in data:
            suboptimal.setdefault(row[0].lower(), row)
    matrices = [
        np.full((len(rows), len(organisms)), None, dtype=object) for _ in range(8)
    ]
    consensus_mrna = [None] * len(rows)
    consensus_srna = [None] * len(rows)
    consensus_names = [None] * len(rows)
    (directory / "evo_alignments").mkdir(exist_ok=True)
    for index, row in enumerate(rows):
        entries = [(j, cell) for j, cell in enumerate(row[2:end]) if cell]
        if len(entries) <= 1:
            continue
        targets = _fasta(directory / f"target_alignments/{index + 1}.aln")
        tabs = []
        subtabs = []
        for j, cell in entries:
            tag = cell.split("(")[0]
            gene = cell.rsplit("(", 1)[-1].split("|")[0].replace("/", "")

            def make(source):
                return [
                    *[
                        int(float(source[k])) if k in (8, 9, 10, 11) else source[k]
                        for k in COORDS
                    ],
                    gene + "_" + tag,
                    tag,
                    float(source[35]),
                    organisms[j],
                ]

            tabs.append(make(optimal[tag]))
            if tag in suboptimal:
                subtabs.append(make(suboptimal[tag]))
        mapped = []
        submapped = []
        for source, destination in ((tabs, mapped), (subtabs, submapped)):
            for tab in source:
                mapped_tab = list(tab)
                tar_map = np.flatnonzero(np.array(list(targets[tab[9]])) != "-") + 1
                srna_map = np.flatnonzero(np.array(list(srna[tab[11]])) != "-") + 1
                for col, positions in (
                    (0, tar_map),
                    (1, tar_map),
                    (4, srna_map),
                    (5, srna_map),
                ):
                    mapped_tab[col] = int(positions[int(tab[col]) - 1])
                destination.append(mapped_tab)
        from .numerics import r_sum

        signal_mrna = np.zeros(
            (len(next(iter(targets.values()))), len(mapped) + len(submapped))
        )
        signal_srna = np.zeros(
            (len(next(iter(srna.values()))), len(mapped) + len(submapped))
        )
        for column, tab in enumerate(mapped + submapped):
            signal_mrna[tab[0] - 1 : tab[1], column] = 1 - tab[10]
            signal_srna[tab[4] - 1 : tab[5], column] = 1 - tab[10]
        signal_mrna = np.array([r_sum(row) for row in signal_mrna])
        signal_srna = np.array([r_sum(row) for row in signal_srna])
        coordinate_strings = {
            k: format_matrix_column([tab[k] for tab in tabs]) for k in (0, 1, 4, 5)
        }
        opt_strings = format_matrix_column([tab[10] for tab in tabs])
        sub_strings = format_matrix_column([tab[10] for tab in subtabs])
        mrna_peak = peak_find(signal_mrna)
        srna_peaks = peak_find_two(signal_srna)
        name = f"{index + 1}_{tabs[0][8]}"
        consensus_names[index] = name
        consensus_mrna[index] = np.array(mrna_peak, dtype=np.int32)
        consensus_srna[index] = [np.array(peak, dtype=np.int32) for peak in srna_peaks]
        outdir = directory / "evo_alignments" / name
        outdir.mkdir(exist_ok=True)

        def label(tab):
            j = organisms.index(tab[11])
            return f"{labels[j]}_{tab[11]}_{tab[9]}"

        _write_fasta(
            outdir / f"{name}_mRNA.fasta",
            [(label(tab), targets[tab[9]]) for tab in tabs],
        )
        _write_fasta(
            outdir / f"{name}_sRNA.fasta", [(label(tab), srna[tab[11]]) for tab in tabs]
        )
        for feature, startcol, endcol, signal, peaks in [
            ("mRNA", 0, 1, signal_mrna, [mrna_peak]),
            ("sRNA", 4, 5, signal_srna, srna_peaks),
        ]:
            annotation = list(COLORS)
            for group in (subtabs, tabs):
                for bucket in range(6):
                    for tab in group:
                        category = int(
                            np.searchsorted(
                                [0.001, 0.01, 0.1, 0.2, 0.3], tab[10], side="left"
                            )
                        )
                        if category == bucket:
                            annotation.append(
                                f"{r_number(tab[10])}\t{label(tab)}\t-1\t{tab[startcol]}\t{tab[endcol]}\t{LABELS[bucket]}"
                            )
            (outdir / f"{name}_{feature}_anno.txt").write_text(
                "\n".join(annotation) + "\n"
            )
            graph = [
                "JALVIEW_ANNOTATION",
                "BAR_GRAPH\tweighted interaction region\t "
                + "|".join(r_number(float(value)) for value in signal),
            ]
            for k, peak in enumerate(peaks):
                if not peak:
                    continue
                text = (
                    "consensus interaction site"
                    if feature == "mRNA"
                    else f"int_site_{k + 1}"
                )
                # Preserve the second sRNA peak's legacy tail-length convention.
                tail = peaks[0][1] if feature == "sRNA" and k == 1 else peak[1]
                site = (
                    "|" * (peak[0] - 1)
                    + "E,"
                    + text
                    + "|".join("E" for _ in range(peak[1] - peak[0]))
                    + "|" * (len(signal) - tail)
                )
                graph.append("NO_GRAPH\t \t " + site)
            (outdir / f"{name}_{feature}_annotation.txt").write_text(
                "\n".join(graph) + "\n"
            )

        def flags(tab, mpeak, speaks):
            return [
                overlap(mpeak, tab[0], tab[1]),
                overlap(speaks[0], tab[4], tab[5]),
                overlap(speaks[1], tab[4], tab[5]),
            ]

        def record(mode, mpeak, speaks):
            output = []
            for tab_number, (tab, aligned) in enumerate(zip(tabs, mapped)):
                j = organisms.index(tab[11])
                booleans = flags(aligned, mpeak, speaks)
                matrices[mode][index, j] = (
                    opt_strings[tab_number]
                    + "|"
                    + "|".join(_boolean(v) for v in booleans)
                )
                sub = next((sub for sub in submapped if sub[8] == tab[8]), None)
                subvalues = [None] * 4
                if sub is not None:
                    subflags = flags(sub, mpeak, speaks)
                    sub_string = sub_strings[
                        next(k for k, t in enumerate(subtabs) if t[8] == tab[8])
                    ]
                    matrices[mode + 1][index, j] = (
                        sub_string + "|" + "|".join(_boolean(v) for v in subflags)
                    )
                    subvalues = [sub_string, *[_boolean(v) for v in subflags]]
                base = [
                    *[
                        coordinate_strings[k][tab_number]
                        if k in coordinate_strings
                        else str(v)
                        if k != 10
                        else opt_strings[tab_number]
                        for k, v in enumerate(tab)
                    ],
                    *[_boolean(v) for v in booleans],
                ]
                output.append(base + subvalues[: 4 if subtabs else 3])
            return output

        result = record(0, mrna_peak, srna_peaks)
        if ooi in [tab[11] for tab in mapped]:
            focus = next(tab for tab in mapped if tab[11] == ooi)
            result = record(2, focus[:2], [focus[4:6], []])
        mapping_header = (
            FIELDS
            + [
                "cons_mRNA",
                "cons_sRNA",
                "cons_sRNA2",
                "p_sub",
                "cons_mRNA_sub",
                "cons_sRNA_sub",
            ]
            + (["cons_sRNA_sub2"] if subtabs else [])
        )
        write_table(
            outdir / f"{name}_mapping_result.txt",
            mapping_header,
            result,
            sep="\t",
            quote=True,
            row_names=[optimal_indices[tab[9]] for tab in tabs],
        )
        for tab in mapped:
            j = organisms.index(tab[11])
            matrices[4][index, j] = f"{tab[0]}|{tab[1]}"
            matrices[6][index, j] = f"{tab[4]}|{tab[5]}"
        for tab in submapped:
            j = organisms.index(tab[11])
            matrices[5][index, j] = f"{tab[0]}|{tab[1]}"
            matrices[7][index, j] = f"{tab[4]}|{tab[5]}"
    save_rdata(
        directory / "conservation_table.Rdata",
        "conservation_table",
        matrices[:4],
        columns=organisms,
    )
    save_rdata(
        directory / "interaction_positions.Rdata",
        "interaction_positions",
        [matrices[i] for i in (4, 5, 6, 7)],
        columns=organisms,
    )
    save_rdata(
        directory / "consensus_positions.Rdata",
        "consensus_both",
        [
            RNamedList(consensus_mrna, consensus_names),
            RNamedList(consensus_srna, consensus_names),
        ],
    )
    return ""
