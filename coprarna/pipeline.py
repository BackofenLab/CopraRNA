"""Complete CopraRNA2 prediction workflow and output organization."""

from dataclasses import dataclass
from pathlib import Path
import logging
import re
import shutil
import zipfile
import tempfile
from . import parameters, organisms, __version__
from .fasta import format_fasta
from .io import ensure_input_is_distinct
from . import homology, plots
from .statistics import r_number

LOG = logging.getLogger("coprarna")


@dataclass
class Config:
    """CopraRNA2 job options; native prediction settings belong in the parameter file.

    Paths supplied by callers are resolved against the job directory. The first
    sRNA record selects the organism of interest unless ``nooi`` is enabled."""

    ntup: int = 200
    ntdown: int = 100
    region: str = "5utr"
    cores: int = 1
    intarna_params: str | None = None
    intarna_params_text: str | None = None
    genome_db: str | None = None
    verbose: bool = False
    websrv: bool = False
    nooi: bool = False
    topcount: int = 200
    root: int = 1
    enrich: int = 0
    noclean: bool = False
    cons: int = 0
    ooifilt: float = 0.0


def _validate(path, config, directory):
    """Validate accessions, sequence content and generated-name collisions.

    Return ordered record IDs without changing the input file. Names that would
    be overwritten by required workflow files are rejected before work begins."""
    if not path.exists():
        raise ValueError("No input FASTA supplied!")
    from Bio import SeqIO

    records = list(SeqIO.parse(path, "fasta"))
    headers = [record.id for record in records]
    if len(headers) < 3:
        raise ValueError("The input file seems to contain less than 3 sequences!")
    if not 0 <= config.ooifilt <= 1:
        raise ValueError("ooifilt needs to be specified between 0 and 1")
    config.cores = max(1, config.cores)
    if config.region not in ("5utr", "3utr", "cds"):
        raise ValueError("region must be one of 5utr, 3utr or cds")
    if config.cons not in (0, 1, 2):
        raise ValueError("cons must be one of 0, 1 or 2")
    if config.cons == 1 and config.nooi:
        raise ValueError("cons 1 can not be combined with nooi")
    if any("-" in str(record.seq) for record in records):
        raise ValueError("Gaps are present in sRNA sequences")
    ensure_input_is_distinct(path, directory)
    for identifier in headers:
        organisms.validate_accession(identifier)
    if len(headers) != len(set(headers)):
        raise ValueError("Duplicate organisms are present in input")
    return headers


def _cleanup(directory, config):
    """Organize generated reports and remove disposable intermediates.

    Guard original FASTA/parameter paths and the saved native-settings snapshot.
    Keep a generated subtree if it contains a protected input. Compressed genomes
    stay in their cache; raw gzip predictions and filter audits are retained."""
    from .runtime import options

    configured = options(directory)
    protected = {
        Path(configured["sRNA FASTA"]).resolve(),
        (directory / "IntaRNA_parameters.txt").resolve(),
    }
    if configured.get("intarna parameter file"):
        protected.add(Path(configured["intarna parameter file"]).resolve())

    def move_generated(path, destination):
        if path.resolve() not in protected and destination.resolve() not in protected:
            shutil.move(str(path), destination)

    patterns = (
        ["ncrna_aligned.fa"]
        + ["cluster.members"]
        + ["aligned_sRNA.fa"]
        + [
            "compatible.distmat.mapped",
            "err.log",
            "*anno*",
            "padj.csv",
            "*pvalues*",
            "ncRNA_*",
            "*.fa.intarna.sorted.csv*",
            "*opt.intarna.csv*",
            "gene_CDS_exception.txt",
            "find_gaps.txt",
            "distmat.out",
            "merged_refseq_ids.txt",
            "CopraRNA2_prep*",
            "fasta_temp_file",
            "fasta_temp_file_out",
            "all.fas*",
            "error.log",
            "formatdb.log",
            "N_chars_in_CDS.txt",
        ]
    )
    if config.enrich:
        patterns += [
            "enrichment.txt",
            "*DAVID*",
            "IntaRNA*",
            "intarna*",
            "org_of_interest_aux_enrichment.txt",
        ]
    if config.websrv:
        patterns += ["*IntaRNA1_ui*", "*top_targets*"]
    for pattern in patterns:
        for path in directory.glob(pattern):
            if path.is_file() and path.resolve() not in protected:
                path.unlink()
    for path in directory.iterdir():
        if (
            re.fullmatch(r".*fa\d+", path.name)
            and path.is_file()
            and path.resolve() not in protected
        ):
            path.unlink()
    for subdir, patterns in [
        ("IntaRNA", ["*.fa.intarna.csv.gz", "*.filter.json"]),
        ("Phylogeny", ["compatible.*", "16s_sequences.*"]),
        ("FASTA", ["*.fa"]),
        ("Regions_plots", ["*regions*"] + (["thumbnail_*"] if config.websrv else [])),
    ]:
        (directory / subdir).mkdir(exist_ok=True)
        for pattern in patterns:
            for path in directory.glob(pattern):
                if path.is_file() and path.resolve() not in protected:
                    move_generated(path, directory / subdir / path.name)
    if config.enrich:
        (directory / "Enrichment").mkdir(exist_ok=True)
        for pattern in (
            "copra_heatmap.html",
            "copraRNA.json",
            "enriched_heatmap_big*",
            "termClusterReport.txt",
            "aux_table.csv",
        ):
            for path in directory.glob(pattern):
                if path.resolve() not in protected:
                    move_generated(path, directory / "Enrichment" / path.name)
    (directory / "all_predictions").mkdir(exist_ok=True)
    for path in directory.glob("*.csv"):
        if path.resolve() not in protected and path.name not in (
            "CopraRNA_result.csv",
            "CopraRNA_result_all.csv",
            "coprarna_websrv_table.csv",
        ):
            move_generated(path, directory / "all_predictions" / path.name)
    (directory / "Rdata").mkdir(exist_ok=True)
    for path in directory.glob("*.Rdata"):
        if path.resolve() not in protected:
            move_generated(path, directory / "Rdata" / path.name)
    for name in ("target_alignments",):
        if not any(directory / name in path.parents for path in protected):
            shutil.rmtree(directory / name, ignore_errors=True)
    evo = directory / "evo_alignments"
    if evo.exists():
        with zipfile.ZipFile(
            directory / "evo_alignments.zip", "w", zipfile.ZIP_DEFLATED
        ) as archive:
            for path in evo.rglob("*"):
                if path.is_file():
                    archive.write(path, path.relative_to(directory))
        if not any(evo in path.parents for path in protected):
            shutil.rmtree(evo)
    warning = directory / "weights.warning"
    if warning.exists() and warning.stat().st_size == 0:
        warning.unlink()


def run_prediction(srna: Path, directory: Path, config: Config | None = None) -> str:
    """Run one complete prediction in an isolated working directory.

    Return empty stdout; scientific results are files and failures propagate as
    exceptions. Preserve the supplied FASTA and parameter file through cleanup.
    """
    return _run_prediction(srna, directory, config or Config())


def _run_prediction(srna, directory, config):
    """Prepare an immutable input copy, database path and scoped native parameters."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    srna = Path(srna)
    srna = srna if srna.is_absolute() else directory / srna
    headers = _validate(srna, config, directory)
    if config.intarna_params is not None and config.intarna_params_text is not None:
        raise ValueError("Supply an IntaRNA parameter file or inline text, not both")
    (directory / "FASTA").mkdir(exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".fa",
        prefix="coprarna_input_",
        dir=directory / "FASTA",
        delete=False,
    ) as stream:
        stream.write(format_fasta(srna))
        working = Path(stream.name)
    database = (
        Path(config.genome_db).expanduser()
        if config.genome_db
        else directory / "genomes"
    )
    if not database.is_absolute():
        database = directory / database
    with parameters.parameter_file(
        directory, config.intarna_params, config.intarna_params_text
    ) as parameter_path:
        return _execute(
            directory,
            config,
            headers,
            srna,
            working,
            database.resolve(),
            parameter_path,
        )


def _execute(directory, config, headers, srna, working, database, parameter_path):
    """Write named options, run homology/prediction, plot results and organize files."""
    version = parameters.check_version(directory)
    labels = [
        ("sRNA FASTA", str(srna.resolve())),
        ("working sRNA FASTA", str(working)),
        ("genome database", str(database)),
        ("intarna parameter file", str(parameter_path) if parameter_path else ""),
        ("IntaRNA version", version),
        ("nt upstream", config.ntup),
        ("nt downstream", config.ntdown),
        ("region", config.region),
        ("core count", config.cores),
        ("verbose", int(config.verbose)),
        ("websrv", int(config.websrv)),
        ("nooi", int(config.nooi)),
        ("top count", config.topcount),
        ("root", config.root),
        ("enrich", config.enrich),
        ("noclean", int(config.noclean)),
        ("cons", config.cons),
        ("ooifilt", config.ooifilt),
        ("version", "CopraRNA " + __version__),
    ]
    (directory / "CopraRNA_option_file.txt").write_text(
        "".join(
            f"{key}:{r_number(value) if isinstance(value, float) else value}\n"
            for key, value in labels
        )
    )
    (directory / "err.log").touch()
    homology.run_homology(
        working,
        config.ntup,
        config.ntdown,
        config.region,
        headers,
        directory,
    )
    warnings = []
    for line in (directory / "zscore.weight").read_text().splitlines():
        identifier, weight = line.split(";")
        if float(weight) > 0.5:
            warnings.append(
                f"{identifier.upper()} may be overweighted. It has weight {float(weight):.2f}. You should consider checking the 16S rDNA tree. We suggest removal of outliers from your input and restarting. "
            )
    (directory / "weights.warning").write_text("".join(warnings))
    if (directory / "err.log").stat().st_size:
        raise RuntimeError("CopraRNA failed. Check err.log for details.")
    try:
        plots.plot_interaction_regions(
            directory / "CopraRNA_result_all.csv", config.topcount, directory
        )
    except ValueError as exc:
        LOG.warning("Interaction-region plotting failed: %s", exc)
    if not config.noclean:
        _cleanup(directory, config)
    return ""
