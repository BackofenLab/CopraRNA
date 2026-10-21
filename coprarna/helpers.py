"""Auditable dispatch for the legacy components that have been ported."""

from dataclasses import dataclass
from pathlib import Path

from .compatibility import legacy_enabled

from . import (
    energy,
    combination,
    annotation,
    exports,
    fasta,
    genome,
    phylogeny,
    resources,
    statistics,
    tables,
)


@dataclass(frozen=True)
class Helper:
    legacy: str
    target: str
    arguments: tuple[str, ...]
    variadic: bool = False


HELPERS = {
    "format-fasta": Helper(
        "format_fasta.pl", "coprarna.fasta.format_fasta", ("fasta",)
    ),
    "fasta-gene": Helper(
        "fasta2genefile.pl", "coprarna.fasta.fasta_to_genefile", ("fasta",)
    ),
    "genome-regions": Helper(
        "parse_region_from_genome.pl",
        "coprarna.genome.parse_regions",
        ("genbank", "ntup", "ntdown", "region"),
    ),
    "check-features": Helper(
        "check_for_gene_CDS_features.pl", "coprarna.genome.check_gene_cds", ("genbank",)
    ),
    "extract-proteins": Helper(
        "get_CDS_from_gbk.pl", "coprarna.genome.extract_proteins", ("genbank",)
    ),
    "extract-16s": Helper(
        "parse_16s_from_gbk.pl",
        "coprarna.genome.extract_16s",
        ("genbank-groups",),
        True,
    ),
    "blast-homology": Helper(
        "blast2homfile.pl", "coprarna.tables.blast_to_homology", ("blast-table",)
    ),
    "disentangle": Helper(
        "disentangle_sorted_intarna_CSV.pl", "coprarna.tables.disentangle", ()
    ),
    "cluster-predictions": Helper(
        "cluster_intarna_csv.pl", "coprarna.tables.cluster_predictions", ()
    ),
    "sampled-count": Helper(
        "get_amount_sampled_values_and_add_to_table.pl",
        "coprarna.tables.add_sampled_count",
        ("result-table", "copra2"),
    ),
    "additional-homologs": Helper(
        "parse_homologs_from_domclust_table.pl",
        "coprarna.tables.additional_homologs",
        ("result-table", "cluster-table"),
    ),
    "enrichment-report": Helper(
        "make_enrichment_table_from_py_output.pl",
        "coprarna.tables.enrichment_report",
        ("service-report",),
    ),
    "transform-distance": Helper(
        "transform_distmat.pl",
        "coprarna.phylogeny.transform_distance",
        ("distance-table",),
    ),
    "compute-weights": Helper(
        "compute_weights.pl",
        "coprarna.phylogeny.compute_weights",
        ("tree", "neighbor-output"),
    ),
    "annotate-intarna": Helper(
        "add_GI_genename_annotation_intarna.pl",
        "coprarna.annotation.annotate_intarna",
        (),
    ),
    "annotate-raw": Helper(
        "annotate_raw_output.pl",
        "coprarna.annotation.annotate_raw",
        ("result-table", "cluster-table", "genbank-groups"),
        True,
    ),
    "adjust-pvalues": Helper(
        "calc_padj.R", "coprarna.statistics.adjust_pvalues", ("pvalue-table",)
    ),
    "prepare-copra2": Helper(
        "prep_cop2.R", "coprarna.statistics.prepare_copra2", ("clustered-table",)
    ),
    "filter-clusters": Helper(
        "remove_clusters_under_percantage.R", "coprarna.statistics.filter_clusters", ()
    ),
    "heatmap-json": Helper(
        "make_heatmap_json.pl", "coprarna.exports.heatmap_json", ("enrichment-table",)
    ),
    "specific-targets": Helper(
        "find_single_specific_targets_in_termCluster.pl",
        "coprarna.exports.specific_targets",
        (),
    ),
    "archive-readme": Helper(
        "print_archive_README.pl", "coprarna.resources.archive_readme", ()
    ),
}

HELPERS.update(
    {
        "fit-energy": Helper(
            "add_pval_to_csv_evdfit.R", "coprarna.energy.add_energy_pvalues", ()
        ),
        "join-copra1": Helper(
            "join_pvals_coprarna1.R",
            "coprarna.combination.join_copra1",
            ("clustered-table",),
        ),
        "join-copra2": Helper(
            "join_pvals_coprarna2.R",
            "coprarna.combination.join_copra2",
            ("organism",),
            True,
        ),
    }
)


HELPERS.update(
    {
        "predict": Helper(
            "CopraRNA2.pl", "coprarna.pipeline.run_prediction", ("srna",), False
        ),
        "homology": Helper(
            "homology_intaRNA.pl",
            "coprarna.homology.run_homology",
            ("srna", "up", "down", "region", "accession"),
            True,
        ),
        "prepare-intarna": Helper(
            "prepare_intarna_out.pl",
            "coprarna.intarna.prepare_predictions",
            ("srna", "up", "down", "region", "genbank-group"),
            True,
        ),
        "combine-clusters": Helper(
            "combine_clusters.pl",
            "coprarna.combination.combine_clusters",
            ("organism-count",),
            False,
        ),
        "download-genbank": Helper(
            "get_refseq_from_refid.pl",
            "coprarna.ncbi.download_genbank",
            ("-acc", "accession", "-g", "output"),
            False,
        ),
        "align-targets": Helper(
            "parallelize_target_alignments.pl",
            "coprarna.alignments.align_targets",
            ("prediction-table",),
            False,
        ),
        "conservation": Helper(
            "copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R",
            "coprarna.conservation.compute_conservation",
            ("organism",),
            False,
        ),
        "plot-regions": Helper(
            "script_R_plots_8.R",
            "coprarna.plots.plot_interaction_regions",
            ("prediction-table", "count"),
            False,
        ),
        "plot-evolution": Helper(
            "evo_heatmap.R", "coprarna.plots.plot_evolutionary_heatmap", (), True
        ),
        "functional-enrichment": Helper(
            "extract_functional_enriched.R",
            "coprarna.enrichment.extract_enriched_terms",
            ("predictions", "report", "output"),
            False,
        ),
        "david-copra": Helper(
            "DAVIDWebService_CopraRNA.py",
            "coprarna.david.coprarna_report",
            ("predictions", "count"),
            False,
        ),
        "david-intarna": Helper(
            "DAVIDWebService_IntaRNA_chartReport.py",
            "coprarna.david.intarna_report",
            ("predictions", "count"),
            False,
        ),
        "webserver-output": Helper(
            "prepare_output_for_websrv_new.pl",
            "coprarna.webserver.prepare_webserver_output",
            ("predictions", "interactions"),
            False,
        ),
        "rasterize": Helper(
            "rasterize.js",
            "coprarna.webserver.rasterize_heatmap",
            ("address", "output"),
            True,
        ),
        "organism-mapping": Helper(
            "build_kegg2refseq.pl",
            "coprarna.maintenance.build_organism_mapping",
            (),
            True,
        ),
        "availability-header": Helper(
            "add_date_omit_incompatible.sh",
            "coprarna.maintenance.write_availability_header",
            ("header-file",),
            False,
        ),
    }
)


def resolve(name: str) -> str:
    if name in HELPERS:
        return name
    for key, helper in HELPERS.items():
        if name == helper.legacy:
            return key
    raise ValueError(
        f'Helper {name!r} has not been migrated; use "coprarna status" for the inventory'
    )


def run_helper(
    name: str, args: list[str], directory: Path, *, seed=None, legacy=False
) -> str:
    from .compatibility import legacy_behavior

    with legacy_behavior(legacy):
        return _run_helper(name, args, directory, seed=seed)


def _run_helper(name: str, args: list[str], directory: Path, *, seed=None) -> str:
    name = resolve(name)
    spec = HELPERS[name]
    expected = len(spec.arguments)
    if len(args) < expected or (not spec.variadic and len(args) != expected):
        raise ValueError(
            f"{name} expects {'at least ' if spec.variadic else ''}{expected} arguments: {' '.join(spec.arguments)}"
        )
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"Working directory does not exist: {directory}")

    def path(index):
        return directory / args[index]

    if name == "format-fasta":
        return fasta.format_fasta(path(0))
    if name == "fasta-gene":
        return fasta.fasta_to_genefile(path(0))
    if name == "genome-regions":
        return genome.parse_regions(path(0), int(args[1]), int(args[2]), args[3])
    if name == "check-features":
        return genome.check_gene_cds(path(0), args[0])
    if name == "extract-proteins":
        return genome.extract_proteins(path(0), resources.organism_mapping(), args[0])
    if name == "extract-16s":
        return genome.extract_16s(args, directory)
    if name == "blast-homology":
        return tables.blast_to_homology(path(0))
    if name == "disentangle":
        return tables.disentangle(directory)
    if name == "cluster-predictions":
        return tables.cluster_predictions(directory)
    if name == "sampled-count":
        return tables.add_sampled_count(path(0), bool(int(args[1])))
    if name == "additional-homologs":
        return tables.additional_homologs(path(0), path(1))
    if name == "enrichment-report":
        return tables.enrichment_report(path(0))
    if name == "transform-distance":
        return phylogeny.transform_distance(path(0))
    if name == "compute-weights":
        return phylogeny.compute_weights(path(0), path(1))
    if name == "annotate-intarna":
        return annotation.annotate_intarna(directory)
    if name == "annotate-raw":
        return annotation.annotate_raw(path(0), path(1), args[2:], directory)
    if name == "adjust-pvalues":
        return statistics.adjust_pvalues(path(0), directory)
    if name == "prepare-copra2":
        return statistics.prepare_copra2(path(0), directory)
    if name == "filter-clusters":
        return statistics.filter_clusters(directory)
    if name == "heatmap-json":
        return exports.heatmap_json(path(0))
    if name == "specific-targets":
        return exports.specific_targets(directory)
    if name == "archive-readme":
        return resources.archive_readme()
    if name == "fit-energy":
        return energy.add_energy_pvalues(directory)
    if name == "join-copra1":
        return combination.join_copra1(path(0), directory, seed)
    if name == "join-copra2":
        return combination.join_copra2(args[0], directory, *args[1:])
    from . import (
        pipeline,
        homology,
        intarna,
        ncbi,
        alignments,
        conservation,
        plots,
        enrichment,
        david,
        webserver,
        maintenance,
    )

    if name == "predict":
        return pipeline.run_prediction(
            path(0),
            directory,
            pipeline.Config(
                seed=seed,
                legacy=legacy_enabled(),
            ),
        )
    if name == "homology":
        return homology.run_homology(
            path(0), int(args[1]), int(args[2]), args[3], args[4:], directory, seed
        )
    if name == "prepare-intarna":
        return intarna.prepare_predictions(
            path(0), int(args[1]), int(args[2]), args[3], args[4:], directory
        )
    if name == "combine-clusters":
        return combination.combine_clusters(int(args[0]), directory, seed)
    if name == "download-genbank":
        if args[0] != "-acc" or args[2] not in ("-g", "--genome"):
            raise ValueError("expected -acc accession -g output")
        return ncbi.download_genbank(args[1], path(3))
    if name == "align-targets":
        return alignments.align_targets(path(0), directory)
    if name == "conservation":
        return conservation.compute_conservation(args[0], directory)
    if name == "plot-regions":
        return plots.plot_interaction_regions(path(0), int(args[1]), directory)
    if name == "plot-evolution":
        return plots.plot_evolutionary_heatmap(
            directory, **dict(arg.split("=", 1) for arg in args)
        )
    if name == "functional-enrichment":
        return enrichment.extract_enriched_terms(path(0), path(1), path(2))
    if name == "david-copra":
        return david.coprarna_report(path(0), int(args[1]))
    if name == "david-intarna":
        return david.intarna_report(path(0), int(args[1]))
    if name == "webserver-output":
        return webserver.prepare_webserver_output(path(0), path(1), directory)
    if name == "rasterize":
        return webserver.rasterize_heatmap(
            args[0],
            path(1),
            args[2] if len(args) > 2 else None,
            float(args[3]) if len(args) > 3 else 1.0,
        )
    if name == "organism-mapping":
        return maintenance.build_organism_mapping(
            directory, path(0) if args else None, seed
        )
    if name == "availability-header":
        return maintenance.write_availability_header(directory, path(0).read_text())
    raise AssertionError(f"Undispatched registered helper: {name}")
