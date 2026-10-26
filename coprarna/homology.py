"""Genome retrieval, homology clustering and comparative prediction workflow."""

from pathlib import Path
import shutil
import json
from Bio import SeqIO

from . import (
    alignments,
    annotation,
    combination,
    conservation,
    fasta,
    genome,
    intarna,
    parameters,
    ncbi,
    resources,
    organisms,
    statistics,
    tables,
)
from .runtime import command, options


def protein_code(index):
    """Return a job-local four-letter DomClust species ID for a zero-based index.

    Native species-name buffers cannot hold full versioned assembly IDs. Four
    lowercase letters keep the validated interface while avoiding a global code
    catalog. Reject negative indices or more than 26**4 organisms; write no files.
    """
    if not 0 <= index < 26**4:
        raise ValueError("DomClust species index must be between 0 and 26**4 - 1")
    return "".join(chr(97 + index // 26**power % 26) for power in (3, 2, 1, 0))


def domclust_table(members, gene_table):
    """Convert supported mode-9 membership output to the removed mode-5 table."""
    species = list(
        dict.fromkeys(
            line.split()[0] for line in gene_table.splitlines() if line.strip()
        )
    )
    entries = []
    domains = {}
    for line in members.splitlines():
        fields = line.split()
        if len(fields) < 4 or ":" not in fields[1]:
            continue
        cluster, identifier = fields[:2]
        domain = fields[2] if len(fields) >= 5 else "1"
        entries.append((cluster, identifier, domain))
        domains.setdefault(identifier, set()).add(domain)
    clusters = {}
    for cluster, identifier, domain in entries:
        organism = identifier.split(":", 1)[0]
        tag = identifier + (f"({domain})" if len(domains[identifier]) > 1 else "")
        clusters.setdefault(cluster, {}).setdefault(organism, []).append(tag)
    if not clusters:
        raise ValueError("DomClust returned no homology clusters")
    return (
        "#id\t"
        + "\t".join(species)
        + "\n"
        + "".join(
            str(i)
            + "\t"
            + "\t".join(" ".join(cells.get(org, [])) for org in species)
            + "\n"
            for i, cells in enumerate(clusters.values(), 1)
        )
    )


def run_domclust(directory):
    """Build cluster.tab from supported DomClust memberships and protein metadata.

    Put flags before positional files and convert supported mode-9 memberships
    to the tabular cluster layout consumed by the Python scientific stages.
    """
    command(
        [
            "domclust",
            "-HO",
            "-S",
            "-c60",
            "-p0.5",
            "-V0.6",
            "-C80",
            "-o9",
            "all.fas.hom",
            "all.fas.gene",
        ],
        directory,
        output="cluster.members",
    )
    table = domclust_table(
        (directory / "cluster.members").read_text(),
        (directory / "all.fas.gene").read_text(),
    )
    (directory / "cluster.tab").write_text(table)


def _finalize_result_table(path):
    """Publish a result CSV with empty missing cells and no internal sort column.

    Rewrite ``path`` in place, converting only complete ``NA`` fields so text that
    happens to contain those letters is unchanged. Return the removed sort-column
    values for the optional webserver alignment map. Reject empty or malformed
    input rather than silently publishing a shifted table.
    """
    lines = path.read_text().splitlines()
    fields = [line.split(",") for line in lines]
    if not fields or any(len(row) < 2 for row in fields):
        raise ValueError(f"Malformed result table: {path}")
    internal_sorting = [row[-1] for row in fields]
    published = [
        ["" if cell == "NA" else cell for cell in row[:-1]] for row in fields
    ]
    text = "\n".join(",".join(row) for row in published) + "\n"
    text = (
        text.replace(",Additional.homologs,", ",Additional homologs,")
        .replace(",Amount.sampled", ",Amount sampled")
        .replace("p.value", "p-value")
    )
    path.write_text(text)
    return internal_sorting


def run_homology(
    srna: Path,
    up: int,
    down: int,
    region: str,
    accessions: list[str],
    directory: Path,
) -> str:
    """Coordinate genomes, protein clusters, native predictions and final reports.

    Resolve replicons once and record their compressed paths for annotation. Build
    protein similarities with BLAST+ and homology with DomClust unless cluster.tab
    was supplied. Combine optimal and suboptimal interactions, map conservation,
    select the requested comparative ranking, then add optional web/enrichment data."""
    parameters.check_version(directory)
    opt = options(directory)
    database = opt.get("genome database") or directory / "genomes"
    metadata = organisms.resolve_many(accessions, database)
    metadata = [
        dict(item, protein_code=protein_code(i)) for i, item in enumerate(metadata)
    ]
    (directory / "genome_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    mapping = {m["query"]: m["replicons"] for m in metadata}
    primary = {
        identifier: replicons[0].split(".")[0]
        for identifier, replicons in mapping.items()
    }
    records = list(SeqIO.parse(srna, "fasta"))
    # All comparative stages index sRNAs by the primary chromosome, including
    # when a user supplies a plasmid, a versioned sequence or an assembly ID.
    (directory / "resolved_sRNA.fa").write_text(
        "".join(f">{primary[record.id]}\n{record.seq}\n" for record in records)
    )
    (directory / "ncrna.fa").write_text(
        "".join(f">ncRNA_{primary[record.id]}\n{record.seq}\n" for record in records)
    )
    groups = []
    codes = {}
    for item in metadata:
        group = [
            ncbi.genome_path(replicon, directory, database)
            for replicon in item["replicons"]
        ]
        groups.append(",".join(str(path) for path in group))
        codes.update({path: item["protein_code"] for path in group})
    genomes = sorted(
        {Path(filename) for group in groups for filename in group.split(",")},
        key=lambda path: path.name,
    )
    (directory / "genome_paths.txt").write_text(
        "".join(str(path) + "\n" for path in genomes)
    )
    exceptions = "".join(genome.check_gene_cds(path, path.name) for path in genomes)
    (directory / "gene_CDS_exception.txt").write_text(exceptions)
    if exceptions:
        raise ValueError("gene but no CDS features present: " + exceptions)
    if not (directory / "cluster.tab").exists():
        with (directory / "all.fas").open("a") as stream:
            for path in genomes:
                stream.write(genome.extract_proteins(path, organism=codes[path]))
        records = list(SeqIO.parse(directory / "all.fas", "fasta"))
        # The legacy uniq -d detects only adjacent repeated FASTA headers.
        duplicates = [a.id for a, b in zip(records, records[1:]) if a.id == b.id]
        (directory / "duplicated_CDS.txt").write_text(
            "".join(">" + tag + "\n" for tag in dict.fromkeys(duplicates))
        )
        if duplicates:
            first = {}
            for record in records:
                if record.id in duplicates:
                    first.setdefault(record.id, record)
            records = [r for r in records if r.id not in duplicates] + list(
                first.values()
            )
            (directory / "all.fas").write_text(
                "".join(f">{r.description}\n{r.seq}\n" for r in records)
            )
        if not (directory / "all.fas.blast").exists():
            command(["makeblastdb", "-in", "all.fas", "-dbtype", "prot"], directory)
            command(
                [
                    "blastp",
                    "-query",
                    "all.fas",
                    "-db",
                    "all.fas",
                    "-evalue",
                    ".001",
                    "-searchsp",
                    "1000000000",
                    "-max_target_seqs",
                    "30000",
                    "-num_threads",
                    opt["core count"],
                    "-outfmt",
                    "6",
                    "-out",
                    "all.fas.blast",
                ],
                directory,
            )
        (directory / "all.fas.hom").write_text(
            tables.blast_to_homology(directory / "all.fas.blast")
        )
        fasta.fasta_to_genefile(directory / "all.fas")
        run_domclust(directory)
    if not (directory / "16s_sequences.fa").exists():
        (directory / "16s_sequences.fa").write_text(
            genome.extract_16s(groups, directory)
        )
    rrna = list(SeqIO.parse(directory / "16s_sequences.fa", "fasta"))
    if len(rrna) != len(accessions) or any("N" in str(record.seq) for record in rrna):
        raise ValueError(
            "wrong number of sequences or N characters in 16s_sequences.fa"
        )
    (directory / "merged_refseq_ids.txt").write_text(" ".join(groups) + "\n")
    intarna.prepare_predictions(
        directory / "ncrna.fa", up, down, region, groups, directory
    )
    combination.combine_clusters(len(accessions), directory)
    prefix, raw, clusters = (
        "CopraRNA2_prep",
        "CopraRNA2_prep_sorted.csv",
        "opt_tags.clustered",
    )
    name = prefix + "_anno.csv"
    (directory / name).write_text(
        annotation.annotate_raw(
            directory / raw, directory / clusters, groups, directory
        )
    )
    name2 = prefix + "_anno_addhomologs.csv"
    (directory / name2).write_text(
        tables.additional_homologs(directory / name, directory / "cluster.tab")
    )
    lines = (directory / name2).read_text().splitlines()
    pvalues = prefix.replace("_prep", "") + "_pvalues.txt"
    (directory / pvalues).write_text(
        "\n".join(line.split(",")[0] for line in lines) + "\n"
    )
    statistics.adjust_pvalues(directory / pvalues, directory)
    adjusted = (directory / "padj.csv").read_text().splitlines()
    name3 = prefix + "_anno_addhomologs_padj.csv"
    (directory / name3).write_text(
        "".join(a + "," + line + "\n" for a, line in zip(adjusted, lines))
    )
    result = prefix + "_anno_addhomologs_padj_amountsamp.csv"
    (directory / result).write_text(tables.add_sampled_count(directory / name3))
    ooi = mapping[accessions[0]][0].split(".")[0]
    alignments.align_targets(
        directory / "CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv", directory
    )
    conservation.compute_conservation(ooi, directory)
    combination.join_copra2(ooi, directory, "ooi_consensus", "overall_consensus")
    suffixes = [
        "_ooi",
        "_balanced",
        "_balanced_consensus",
        "_ooi_consensus",
        "_ooi_ooiconsensus",
    ]
    for suffix in suffixes:
        base = "CopraRNA2"
        lines = (
            (directory / (base + "_final_all" + suffix + ".csv"))
            .read_text()
            .splitlines(keepends=True)
        )
        (directory / (base + "_final" + suffix + ".csv")).write_text(
            "".join(lines[: int(opt["top count"]) + 1])
        )
    if int(opt["nooi"]):
        chosen = "CopraRNA2_final_balanced" + (
            "_consensus" if int(opt["cons"]) == 2 else ""
        )
    else:
        chosen = "CopraRNA2_final_ooi" + (
            "_ooiconsensus"
            if int(opt["cons"]) == 1
            else "_consensus"
            if int(opt["cons"]) == 2
            else ""
        )
    all_name = chosen.replace("_final", "_final_all")
    for source, target in [
        (chosen + ".csv", "CopraRNA_result.csv"),
        (all_name + ".csv", "CopraRNA_result_all.csv"),
    ]:
        shutil.copyfile(directory / source, directory / target)
    if float(opt["ooifilt"]):
        lines = (
            (directory / "CopraRNA_result_all.csv")
            .read_text()
            .splitlines(keepends=True)
        )

        def keep(line):
            cell = line.split(",")[2]
            return cell and float(cell.split("|")[2]) <= float(opt["ooifilt"])

        lines = (
            lines[:1]
            + [line for line in lines[1:] if keep(line)]
            + [line for line in lines[1:] if not keep(line)]
        )
        (directory / "CopraRNA_result_all_filt.csv").write_text("".join(lines))
        (directory / "CopraRNA_result_all.csv").write_text("".join(lines))
        (directory / "CopraRNA_result.csv").write_text(
            "".join(lines[: int(opt["top count"]) + 1])
        )
    from .plots import plot_evolutionary_heatmap

    plot_evolutionary_heatmap(directory)
    for filename in ("CopraRNA_result.csv", "CopraRNA_result_all.csv"):
        path = directory / filename
        internal_sorting = _finalize_result_table(path)
        if filename == "CopraRNA_result.csv" and int(opt["websrv"]):
            (directory / "CopraRNA_result.map_evo_align").write_text(
                "\n".join(internal_sorting) + "\n"
            )
    if int(opt["websrv"]):
        from .webserver import export_prediction

        export_prediction(directory, ooi, up, down)
    (directory / "README.txt").write_text(resources.archive_readme())
    if int(opt["enrich"]):
        from .enrichment import enrich_prediction

        enrich_prediction(directory, int(opt["enrich"]))
    return ""
