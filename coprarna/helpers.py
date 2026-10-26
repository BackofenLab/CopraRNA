"""Named Python helpers for individual scientific workflow stages."""

from dataclasses import dataclass
from pathlib import Path


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
    """Describe a Python helper and its positional CLI argument contract.

    Variadic helpers accept trailing arguments after the required argument names.
    """

    target: str
    arguments: tuple[str, ...]
    variadic: bool = False


HELPERS = {
    "format-fasta": Helper("coprarna.fasta.format_fasta", ("fasta",)),
    "fasta-gene": Helper("coprarna.fasta.fasta_to_genefile", ("fasta",)),
    "genome-regions": Helper(
        "coprarna.genome.parse_regions",
        ("genbank", "ntup", "ntdown", "region"),
    ),
    "check-features": Helper("coprarna.genome.check_gene_cds", ("genbank",)),
    "extract-proteins": Helper("coprarna.genome.extract_proteins", ("genbank",)),
    "extract-16s": Helper(
        "coprarna.genome.extract_16s",
        ("genbank-groups",),
        True,
    ),
    "blast-homology": Helper("coprarna.tables.blast_to_homology", ("blast-table",)),
    "disentangle": Helper("coprarna.tables.disentangle", ()),
    "cluster-predictions": Helper("coprarna.tables.cluster_predictions", ()),
    "sampled-count": Helper(
        "coprarna.tables.add_sampled_count",
        ("result-table",),
    ),
    "additional-homologs": Helper(
        "coprarna.tables.additional_homologs",
        ("result-table", "cluster-table"),
    ),
    "enrichment-report": Helper(
        "coprarna.tables.enrichment_report",
        ("service-report",),
    ),
    "transform-distance": Helper(
        "coprarna.phylogeny.transform_distance",
        ("distance-table",),
    ),
    "compute-weights": Helper(
        "coprarna.phylogeny.compute_weights",
        ("tree", "neighbor-output"),
    ),
    "annotate-intarna": Helper(
        "coprarna.annotation.annotate_intarna",
        (),
    ),
    "annotate-raw": Helper(
        "coprarna.annotation.annotate_raw",
        ("result-table", "cluster-table", "genbank-groups"),
        True,
    ),
    "adjust-pvalues": Helper("coprarna.statistics.adjust_pvalues", ("pvalue-table",)),
    "prepare-copra2": Helper(
        "coprarna.statistics.prepare_copra2", ("clustered-table",)
    ),
    "heatmap-json": Helper("coprarna.exports.heatmap_json", ("enrichment-table",)),
    "specific-targets": Helper(
        "coprarna.exports.specific_targets",
        (),
    ),
    "archive-readme": Helper("coprarna.resources.archive_readme", ()),
}

HELPERS.update(
    {
        "fit-energy": Helper("coprarna.energy.add_energy_pvalues", ()),
        "join-copra2": Helper(
            "coprarna.combination.join_copra2",
            ("organism",),
            True,
        ),
    }
)


HELPERS.update(
    {
        "predict": Helper("coprarna.pipeline.run_prediction", ("srna",), False),
        "homology": Helper(
            "coprarna.homology.run_homology",
            ("srna", "up", "down", "region", "accession"),
            True,
        ),
        "prepare-intarna": Helper(
            "coprarna.intarna.prepare_predictions",
            ("srna", "up", "down", "region", "genbank-group"),
            True,
        ),
        "combine-clusters": Helper(
            "coprarna.combination.combine_clusters",
            ("organism-count",),
            False,
        ),
        "download-genbank": Helper(
            "coprarna.ncbi.download_genbank",
            ("-acc", "accession", "-g", "output"),
            False,
        ),
        "align-targets": Helper(
            "coprarna.alignments.align_targets",
            ("prediction-table",),
            False,
        ),
        "conservation": Helper(
            "coprarna.conservation.compute_conservation",
            ("organism",),
            False,
        ),
        "plot-regions": Helper(
            "coprarna.plots.plot_interaction_regions",
            ("prediction-table", "count"),
            False,
        ),
        "plot-evolution": Helper("coprarna.plots.plot_evolutionary_heatmap", (), True),
        "functional-enrichment": Helper(
            "coprarna.enrichment.extract_enriched_terms",
            ("predictions", "report", "output"),
            False,
        ),
        "david-copra": Helper(
            "coprarna.david.coprarna_report",
            ("predictions", "count"),
            False,
        ),
        "david-intarna": Helper(
            "coprarna.david.intarna_report",
            ("predictions", "count"),
            False,
        ),
        "webserver-output": Helper(
            "coprarna.webserver.prepare_webserver_output",
            ("predictions", "interactions"),
            False,
        ),
        "rasterize": Helper(
            "coprarna.webserver.rasterize_heatmap",
            ("address", "output"),
            True,
        ),
    }
)


def resolve(name: str) -> str:
    """Resolve an active Python helper name or raise ValueError."""
    if name in HELPERS:
        return name
    raise ValueError(
        f'Unknown helper {name!r}; use "coprarna helpers" for available names'
    )


def run_helper(name: str, args: list[str], directory: Path) -> str:
    """Run a named helper in an existing directory, returning its text output.

    Validate the positional contract before dispatch. Helper-specific file outputs
    are written in directory; input, native-tool and scientific errors propagate.
    """
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
        return genome.extract_proteins(path(0))
    if name == "extract-16s":
        return genome.extract_16s(args, directory)
    if name == "blast-homology":
        return tables.blast_to_homology(path(0))
    if name == "disentangle":
        return tables.disentangle(directory)
    if name == "cluster-predictions":
        return tables.cluster_predictions(directory)
    if name == "sampled-count":
        return tables.add_sampled_count(path(0))
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
    if name == "heatmap-json":
        return exports.heatmap_json(path(0))
    if name == "specific-targets":
        return exports.specific_targets(directory)
    if name == "archive-readme":
        return resources.archive_readme()
    if name == "fit-energy":
        return energy.add_energy_pvalues(directory)
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
    )

    if name == "predict":
        return pipeline.run_prediction(
            path(0),
            directory,
            pipeline.Config(),
        )
    if name == "homology":
        return homology.run_homology(
            path(0), int(args[1]), int(args[2]), args[3], args[4:], directory
        )
    if name == "prepare-intarna":
        return intarna.prepare_predictions(
            path(0), int(args[1]), int(args[2]), args[3], args[4:], directory
        )
    if name == "combine-clusters":
        return combination.combine_clusters(int(args[0]), directory)
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
    raise AssertionError(f"Undispatched registered helper: {name}")
