"""Streaming access to compressed and plain scientific data files."""

from contextlib import contextmanager
from fnmatch import fnmatchcase
import gzip
from pathlib import Path
from Bio import SeqIO


@contextmanager
def open_text(path, mode="rt"):
    """Open UTF-8 text directly from plain or gzip storage; preserve CSV newlines."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, mode, encoding="utf-8", newline="") as stream:
        yield stream


def read_text(path):
    """Read a small complete text table through the shared compression adapter."""
    with open_text(path) as stream:
        return stream.read()


def write_text(path, text):
    """Write a text table as gzip when its destination has a .gz suffix."""
    with open_text(path, "wt") as stream:
        stream.write(text)


def paths(directory, pattern):
    """Prefer compressed data when both representations of a table exist."""
    files = {p.name: p for p in Path(directory).glob(pattern)}
    for path in Path(directory).glob(pattern + ".gz"):
        files[path.name[:-3]] = path
    return [files[name] for name in sorted(files)]


def accession(path):
    """Extract the RefSeq label from plain/compressed GenBank file names."""
    name = Path(path).name.removesuffix(".gz")
    for suffix in (".gbff", ".gbk", ".gb"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def genbank_records(path):
    """Yield GenBank records from the compressed stream without a plain disk copy."""
    with open_text(path) as stream:
        yield from SeqIO.parse(stream, "genbank")


def genome_paths(directory):
    """Read the job genome manifest, or discover local genomes for standalone helpers."""
    manifest = Path(directory) / "genome_paths.txt"
    if manifest.exists():
        return [Path(line) for line in manifest.read_text().splitlines() if line]
    return sorted(
        {
            p
            for suffix in ("*.gb", "*.gb.gz", "*.gbk.gz", "*.gbff.gz")
            for p in Path(directory).glob(suffix)
        }
    )


# Fixed/native output namespaces. Inputs may use any other basename, including
# input_sRNA.fa and custom.csv; keep this list current when adding new producers.
_JOB_OUTPUTS = (
    "ncrna.fa",
    "ncrna_aligned.fa",
    "aligned_sRNA.fa",
    "ncRNA_*.fa",
    "*_upfromstartpos_*_down_*.fa*",
    "*.filter.json",
    "16s_sequences.*",
    "compatible.*",
    "zscore.weight",
    "weights.warning",
    "genome_paths.txt",
    "CopraRNA_option_file.txt",
    "CopraRNA*.csv",
    "CopraRNA*_pvalues.txt",
    "CopraRNA_top_targets.fa",
    "CopraRNA_result.map_evo_align",
    "IntaRNA_parameters.txt",
    "IntaRNA_web.*",
    "IntaRNA_chartReport*.txt",
    "intarna_websrv_table*.csv",
    "coprarna_*table.csv",
    "all.fas*",
    "cluster.tab",
    "cluster.members",
    "opt_tags.*",
    "subopt_tags.*",
    "padj.csv",
    "fasta_temp_file*",
    "distmat*",
    "err.log",
    "error.log",
    "README.txt",
    "gene_CDS_exception.txt",
    "duplicated_CDS.txt",
    "merged_refseq_ids.txt",
    "conservation_table.Rdata",
    "interaction_positions.Rdata",
    "consensus_positions.Rdata",
    "evo_alignments.zip",
    "mRNA_regions*",
    "sRNA_regions*",
    "thumbnail_*.svg",
    "*_conservation_heatmap.pdf",
    "*_conservation_heatmap.svg",
    "enriched_heatmap_big.*",
    "enrichment.txt",
    "DAVID_enrichment*.txt",
    "termClusterReport.txt",
    "aux_table.csv",
    "org_of_interest_aux_enrichment.txt",
    "copra_heatmap.html",
    "copraRNA.json",
)


def ensure_input_is_distinct(path, directory, *, allow_parameter_snapshot=False):
    """Reject inputs that a required job output would overwrite or misinterpret.

    Check both the supplied path and canonical path, plus existing output symlinks
    that resolve to the input. Check before writing potentially conflicting files.
    A supplied native
    settings snapshot is allowed because it is forwarded without rewriting it.
    Cleanup separately guards destinations when moving reports into subfolders.
    """
    directory = Path(directory).resolve()
    original = Path(path).resolve()
    candidates = {Path(path).absolute(), original}
    candidates.update(
        entry
        for entry in directory.iterdir()
        if entry.is_symlink() and entry.resolve() == original
    )
    for candidate in candidates:
        try:
            relative = candidate.relative_to(directory)
        except ValueError:
            continue
        if len(relative.parts) > 1:
            if relative.parts[0] in {"target_alignments", "evo_alignments"}:
                raise ValueError(
                    f"Input {path} collides with the generated {relative.parts[0]} directory; "
                    "use a separate working directory"
                )
            continue
        if allow_parameter_snapshot and candidate.name == "IntaRNA_parameters.txt":
            continue
        if any(fnmatchcase(candidate.name, pattern) for pattern in _JOB_OUTPUTS):
            raise ValueError(
                f"Input {path} collides with generated output {candidate.name}; "
                "use a separate working directory"
            )
