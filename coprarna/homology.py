"""Genome retrieval, homology clustering and comparative prediction workflow."""

from pathlib import Path
import shutil
import time
from Bio import SeqIO

from . import (
    alignments,
    annotation,
    combination,
    conservation,
    fasta,
    genome,
    intarna,
    ncbi,
    resources,
    statistics,
    tables,
)
from .runtime import command, options


def affiliations():
    result = {}
    for line in resources.organism_mapping().splitlines():
        if "\t" not in line:
            continue
        code, ids = line.split("\t")
        for identifier in ids.split():
            result[identifier] = ids.split()
    return result


def run_domclust(directory):
    from .compatibility import legacy_enabled

    flags = ["-HO", "-S", "-c60", "-p0.5", "-V0.6", "-C80", "-o5"]
    files = ["all.fas.hom", "all.fas.gene"]
    command(
        ["domclust", *(files + flags if legacy_enabled() else flags + files)],
        directory,
        output="cluster.tab",
    )


def run_homology(
    srna: Path,
    up: int,
    down: int,
    region: str,
    accessions: list[str],
    directory: Path,
    seed=None,
) -> str:
    opt = options(directory)
    mapping = affiliations()
    groups = []
    text = srna.read_text().replace(">", ">ncRNA_")
    for identifier in accessions:
        text = text.replace(identifier, mapping[identifier][0])
    (directory / "ncrna.fa").write_text(text)
    for identifier in accessions:
        filenames = []
        for replicon in mapping[identifier]:
            filename = replicon + ".gb"
            filenames.append(filename)
            ncbi.download_genbank(replicon, directory / filename)
        groups.append(",".join(filenames))
    for path in sorted(directory.glob("*gb")):
        lines = path.read_text().splitlines(keepends=True)
        if not any(line.startswith("//") for line in lines[-2:]):
            path.unlink()
            delay = 30.0
            for _ in range(101):
                time.sleep(delay)
                delay *= 1.1
                ncbi.download_genbank(path.name[:-3], path)
                if path.exists() and any(
                    line.startswith("//") for line in path.read_text().splitlines()[-2:]
                ):
                    break
            else:
                raise RuntimeError(
                    "Not all RefSeq *gb files downloaded correctly. Restart your job."
                )
        path.write_text(
            "".join(
                line
                for line in path.read_text().splitlines(keepends=True)
                if not line.startswith("CONTIG")
            )
        )
    exceptions = "".join(
        genome.check_gene_cds(path, path.name) for path in sorted(directory.glob("*gb"))
    )
    (directory / "gene_CDS_exception.txt").write_text(exceptions)
    if exceptions:
        raise ValueError("gene but no CDS features present: " + exceptions)
    if not (directory / "cluster.tab").exists():
        with (directory / "all.fas").open("a") as stream:
            for path in sorted(directory.glob("*gb")):
                stream.write(
                    genome.extract_proteins(
                        path, resources.organism_mapping(), path.name
                    )
                )
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
            command(["formatdb", "-i", "all.fas"], directory)
            command(
                [
                    "blastall",
                    "-a",
                    opt["core count"],
                    "-p",
                    "blastp",
                    "-d",
                    "all.fas",
                    "-e",
                    ".001",
                    "-i",
                    "all.fas",
                    "-Y",
                    "1e9",
                    "-v",
                    "30000",
                    "-b",
                    "30000",
                    "-m",
                    "8",
                    "-o",
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
    combination.combine_clusters(len(accessions), directory, seed)
    cop1 = bool(int(opt["CopraRNA1"]))
    for prefix, raw, clusters in [
        ("CopraRNA2_prep", "CopraRNA2_prep_sorted.csv", "opt_tags.clustered")
    ] + (
        [
            (
                "CopraRNA1",
                "CopraRNA1_with_pvsample_sorted.csv",
                "opt_tags.clustered_rcsize",
            )
        ]
        if cop1
        else []
    ):
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
        (directory / result).write_text(
            tables.add_sampled_count(directory / name3, prefix != "CopraRNA1")
        )
        if prefix == "CopraRNA1":
            (directory / result).rename(directory / "CopraRNA1_final_all.csv")
    ooi = mapping[accessions[0]][0]
    if not cop1:
        alignments.align_targets(
            directory / "CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv", directory
        )
        conservation.compute_conservation(ooi, directory)
        combination.join_copra2(ooi, directory, "ooi_consensus", "overall_consensus")
    suffixes = (
        [""]
        if cop1
        else [
            "_ooi",
            "_balanced",
            "_balanced_consensus",
            "_ooi_consensus",
            "_ooi_ooiconsensus",
        ]
    )
    for suffix in suffixes:
        base = "CopraRNA1" if cop1 else "CopraRNA2"
        lines = (
            (directory / (base + "_final_all" + suffix + ".csv"))
            .read_text()
            .splitlines(keepends=True)
        )
        (directory / (base + "_final" + suffix + ".csv")).write_text(
            "".join(lines[: int(opt["top count"]) + 1])
        )
    if cop1:
        chosen = "CopraRNA1_final"
    elif int(opt["nooi"]):
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
    if not cop1:
        from .plots import plot_evolutionary_heatmap

        plot_evolutionary_heatmap(directory)
        for filename in ("CopraRNA_result.csv", "CopraRNA_result_all.csv"):
            path = directory / filename
            lines = path.read_text().splitlines()
            if filename == "CopraRNA_result.csv" and int(opt["websrv"]):
                (directory / "CopraRNA_result.map_evo_align").write_text(
                    "\n".join(line.split(",")[-1] for line in lines) + "\n"
                )
            text = "\n".join(",".join(line.split(",")[:-1]) for line in lines) + "\n"
            text = (
                text.replace(",Additional.homologs,", ",Additional homologs,")
                .replace(",Amount.sampled", ",Amount sampled")
                .replace("p.value", "p-value")
            )
            path.write_text(text)
    if int(opt["websrv"]):
        from .webserver import export_prediction

        export_prediction(directory, ooi, up, down)
    (directory / "README.txt").write_text(resources.archive_readme())
    if int(opt["enrich"]):
        from .enrichment import enrich_prediction

        enrich_prediction(directory, int(opt["enrich"]))
    return ""
