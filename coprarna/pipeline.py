"""Complete Python command-line workflow, retaining CopraRNA 2.1.3 options."""

from dataclasses import dataclass
from pathlib import Path
import logging
import re
import shutil
import zipfile
from .fasta import format_fasta
from . import homology, plots
from .statistics import r_number

LOG = logging.getLogger("coprarna")


@dataclass
class Config:
    ntup: int = 200
    ntdown: int = 100
    region: str = "5utr"
    rcsize: float = 0.5
    cores: int = 1
    winsize: int = 150
    maxbpdist: int = 100
    cop1: bool = False
    verbose: bool = False
    websrv: bool = False
    nooi: bool = False
    topcount: int = 200
    root: int = 1
    enrich: int = 0
    noclean: bool = False
    cons: int = 0
    ooifilt: float = 0.0
    seed: int | None = None


def _validate(path, config, directory):
    if not path.exists():
        raise ValueError("No input FASTA supplied!")
    text = path.read_text()
    headers = [
        line.replace(">", "").strip() for line in text.splitlines() if ">" in line
    ]
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
    if config.cons and config.cop1:
        raise ValueError("cons must be 0 for cop1 prediction")
    (directory / "find_gaps.txt").write_text(
        "".join(line + "\n" for line in text.splitlines() if "-" in line)
    )
    if "-" in text:
        raise ValueError("Gaps are present in sRNA sequences")
    mapping = homology.affiliations()
    for identifier in headers:
        if not re.search(r"NC_\d{6}|NZ_.*", identifier):
            raise ValueError(f"{identifier} does not match correct RefSeq ID format")
        if identifier not in mapping:
            raise ValueError(f"{identifier} is not present in the availability list")
    if config.maxbpdist > config.winsize:
        raise ValueError("maxbpdist must be <= winsize")
    if (
        config.region in ("5utr", "3utr")
        and config.ntup + config.ntdown < config.winsize
    ):
        raise ValueError("ntup + ntdown must be >= winsize")
    if not 0 <= config.rcsize <= 1:
        raise ValueError("rcsize must be between 0 and 1")
    if len(headers) != len(set(headers)):
        raise ValueError("Duplicate organisms are present in input")
    return headers


def _cleanup(directory, config):
    patterns = (
        ["ncrna_aligned.fa"]
        + ([] if config.cop1 else ["aligned_sRNA.fa"])
        + [
            "*.gb",
            "compatible.distmat.mapped",
            "err.log",
            "*anno*",
            "padj.csv",
            "*pvalues*",
            "ncRNA_*",
            "*.fa.intarna.sorted.csv",
            "*opt.intarna.csv",
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
    if config.cop1:
        patterns += ["*pvsample*", "CopraRNA1_final_all.csv", "CopraRNA1_final.csv"]
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
            if path.is_file():
                path.unlink()
    for path in directory.iterdir():
        if re.fullmatch(r".*fa\d+", path.name) and path.is_file():
            path.unlink()
    for subdir, patterns in [
        ("IntaRNA", ["*.fa.intarna.csv"]),
        ("Phylogeny", ["compatible.*", "16s_sequences.*"]),
        ("FASTA", ["*.fa"]),
        ("Regions_plots", ["*regions*"] + (["thumbnail_*"] if config.websrv else [])),
    ]:
        (directory / subdir).mkdir(exist_ok=True)
        for pattern in patterns:
            for path in directory.glob(pattern):
                if path.is_file():
                    shutil.move(str(path), directory / subdir / path.name)
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
                shutil.move(str(path), directory / "Enrichment" / path.name)
    if not config.cop1:
        (directory / "all_predictions").mkdir(exist_ok=True)
        for path in directory.glob("*.csv"):
            if path.name not in (
                "CopraRNA_result.csv",
                "CopraRNA_result_all.csv",
                "coprarna_websrv_table.csv",
            ):
                shutil.move(str(path), directory / "all_predictions" / path.name)
        (directory / "Rdata").mkdir(exist_ok=True)
        for path in directory.glob("*.Rdata"):
            shutil.move(str(path), directory / "Rdata" / path.name)
        for name in ("target_alignments",):
            shutil.rmtree(directory / name, ignore_errors=True)
    evo = directory / "evo_alignments"
    if evo.exists():
        with zipfile.ZipFile(
            directory / "evo_alignments.zip", "w", zipfile.ZIP_DEFLATED
        ) as archive:
            for path in evo.rglob("*"):
                if path.is_file():
                    archive.write(path, path.relative_to(directory))
        shutil.rmtree(evo)
    warning = directory / "weights.warning"
    if warning.exists() and warning.stat().st_size == 0:
        warning.unlink()


def run_prediction(srna: Path, directory: Path, config: Config | None = None) -> str:
    config = config or Config()
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    srna = Path(srna)
    srna = srna if srna.is_absolute() else directory / srna
    headers = _validate(srna, config, directory)
    labels = [
        ("sRNA FASTA", str(srna) if srna.name != "input_sRNA.fa" else "input_sRNA.fa"),
        ("nt upstream", config.ntup),
        ("nt downstream", config.ntdown),
        ("region", config.region),
        ("relative clustersize", config.rcsize),
        ("core count", config.cores),
        ("win size", config.winsize),
        ("max bp dist", config.maxbpdist),
        ("CopraRNA1", int(config.cop1)),
        ("verbose", int(config.verbose)),
        ("websrv", int(config.websrv)),
        ("nooi", int(config.nooi)),
        ("top count", config.topcount),
        ("root", config.root),
        ("enrich", config.enrich),
        ("noclean", int(config.noclean)),
        ("cons", config.cons),
        ("ooifilt", config.ooifilt),
        ("version", "CopraRNA 2.1.3"),
    ]
    (directory / "CopraRNA_option_file.txt").write_text(
        "".join(
            f"{key}:{r_number(value) if isinstance(value, float) else value}\n"
            for key, value in labels
        )
    )
    formatted = format_fasta(srna)
    (directory / "input_sRNA.fa").write_text(formatted)
    (directory / "err.log").touch()
    homology.run_homology(
        directory / "input_sRNA.fa",
        config.ntup,
        config.ntdown,
        config.region,
        headers,
        directory,
        config.seed,
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
        LOG.warning("Legacy-compatible plotting failure: %s", exc)
    if not config.noclean:
        _cleanup(directory, config)
    return ""
