# CopraRNA 3

**Comparative target prediction for prokaryotic trans-acting small RNAs.**

CopraRNA combines whole-genome IntaRNA predictions across homologous sRNAs from
at least three distinct organisms. Phylogenetic weighting, energy statistics and
optional interaction-site conservation produce a ranked list of candidate targets.
The first FASTA record selects the organism of interest by default.

Version **3.0.0** is a Python 3 implementation of the CopraRNA2 comparative
algorithm. It supports IntaRNA 3, compressed shared genome databases, independent
batches, native parameter files, PDF/SVG figures and Jalview alignments. See
[ChangeLog.md](ChangeLog.md) for the transition from 2.1.3 and the scientific fixes.

- [Installation](#installation)
- [Input and a first prediction](#input-and-a-first-prediction)
- [Genome database and supported accessions](#genome-database-and-supported-accessions)
- [Prediction options](#prediction-options)
- [Outputs and interpretation](#outputs-and-interpretation)
- [Independent batches and helpers](#independent-batches-and-helpers)
- [Upgrading from 2.1](#upgrading-from-21)
- [Development and validation](#development-and-validation)
- [Citation](#citation)

## Installation

### Complete scientific environment

Python **3.11 or newer** is required. The tested complete environment targets
Linux x86-64 and uses Python 3.12. From a checkout of this release:

```bash
git clone https://github.com/BackofenLab/CopraRNA.git
cd CopraRNA
# Select the 3.0 release checkout before installing.
conda env create -f environment.yml
conda activate coprarna-python3
python -m pip install --no-deps .
coprarna --version
coprarna run --help
```

The version command should print `3.0.0`. Until this branch is merged, check out
`codex/issue-4-release-3` to try the release preparation. The environment file is
not a cross-platform lockfile. Existing CopraRNA 2 Bioconda packages and container
images do not install this Python release.

| Dependency | Validated environment / purpose |
| --- | --- |
| IntaRNA / ViennaRNA | 3.4.1 / 2.7.2; `IntaRNAsTar` predictions and native gzip CSV |
| BLAST+ | 2.17.0; `makeblastdb` and `blastp` protein similarities |
| DomClust | 1.2.8; protein homology clustering, membership output mode 9 |
| EMBOSS / embassy-phylip | 6.6.0 / 3.69.660; `distmat` and `fneighbor` |
| MAFFT / Clustal Omega | 7.526 / 1.2.4; target and sRNA alignments |
| Biopython | 1.88; parsing, sequences and trees; 1.84 also tested on Python 3.11 |
| NumPy / SciPy | 2.5.4 / 1.18.1 on Python >=3.12; 2.4.3 / 1.17.1 on Python 3.11 |
| pandas / Matplotlib | 3.0.6 / 3.11.2; scientific data and vector reports |
| rdata / suds-community | 1.1.0 / 1.2.0; scientific matrix serialization and optional DAVID requests |

The Python controller does not require Perl, R, PhantomJS, PostScript or
Ghostscript. Native packages can have their own transitive interpreter
dependencies. Matplotlib is installed from its PyPI wheel to avoid conflicting
ICU requirements between conda rendering packages and IntaRNA.

### Python package and existing tools

For Python helpers or an existing native-tool installation:

```bash
python -m pip install .
export COPRARNA_TOOL_PATH=/opt/scientific-tools/bin
```

`COPRARNA_TOOL_PATH` prepends one or more colon-separated directories to the
subprocess search path. Complete predictions need all the native executables
listed above, including `IntaRNA` and `IntaRNAsTar`. IntaRNA releases older than 3
are rejected. The pinned numerical versions preserve validated scientific fits;
changing those pins requires new validation.

## Input and a first prediction

Prepare an **unaligned nucleotide FASTA** containing at least three homologous
sRNAs from three distinct assemblies. Use an NCBI RefSeq chromosome/plasmid
accession (`NC_...` or `NZ_...`) or an assembly accession (`GCF_...`) as each record
ID. Optional version suffixes are supported. Descriptions after the ID are ignored.
Sequences must not contain alignment gaps (`-`). Select homologous sRNAs and a
useful phylogenetic spread; including nearly identical organisms can reduce the
information contributed by comparative prediction.

For example, the record IDs for E. coli K-12 MG1655, B. subtilis 168 and
P. aeruginosa PAO1 are `NC_000913`, `NC_000964` and `NC_002516`. Supply the actual
homologous sRNA sequences for your study; genome accessions are identifiers, not
sRNA sequences. Put the organism of interest first. A second chromosome or plasmid
from the same assembly is not an independent organism and is rejected.

Use a dedicated job directory. Relative FASTA and parameter-file paths are
resolved against **the job directory**, not the shell's current directory:

```bash
mkdir -p /data/jobs/my-srna
cp /data/my_sRNAs.fa /data/jobs/my-srna/input_sRNA.fa
coprarna run --workdir /data/jobs/my-srna --cores 8 \
  --genome-db /data/coprarna-genomes
```

Or use an absolute input path without copying it:

```bash
coprarna run --workdir /data/jobs/another-srna \
  --srnaseq /data/my_sRNAs.fa --genome-db /data/coprarna-genomes
```

The original FASTA remains byte unchanged, including when named `input_sRNA.fa`.
A normalized working copy goes under `FASTA/`. Inputs that collide with required
job outputs, including symlink aliases and protected cleanup destinations, are
rejected. Choose a new job directory when changing inputs, genome versions or
clustering settings: intermediate homology files are reusable within a job, and
changing a completed job in place is not a full recomputation mechanism.

Runtime and memory depend on genome sizes, organism count, RNA lengths and native
settings. Whole-genome predictions can take hours; start with at least 8 GB RAM
and allow disk space for protein searches, alignments and reports.

## Genome database and supported accessions

### On-demand resolution replaces availability lists

There is no bundled organism allow-list to regenerate. CopraRNA resolves the
requested IDs through [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
and the [NCBI Datasets v2 API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/).
It obtains the current complete RefSeq assembly, organism name and **all** of its
chromosomes and plasmids. Newly available assemblies do not require a CopraRNA
release or a KEGG-to-RefSeq table update.

Supported inputs resolve to one current **Complete Genome** RefSeq assembly with
complete `NC_`/`NZ_` replicons. Genome annotations must provide translated CDSs
with locus tags and an identifiable 16S rRNA. Draft contig/scaffold assemblies,
missing annotations and ambiguous assembly membership produce diagnostics. If a
sequence belongs to multiple current complete assemblies, choose the intended
`GCF_` assembly as the FASTA ID. Explicit old sequence versions are not silently
substituted with a different version in a current assembly.

Inspect the resolution before a long prediction:

```bash
coprarna genomes NC_000913 NC_003062 --genome-db /data/coprarna-genomes
coprarna genomes GCF_000005845.2 --genome-db /data/coprarna-genomes
```

The JSON output includes the query ID, exact assembly version, organism name,
ordered versioned replicons and retrieval time. The longest chromosome labels
an organism's report columns; other chromosomes and plasmids remain included in
prediction. Reports use unversioned accession labels, while cache files and job
metadata record exact sequence versions.

### Prefetch, reuse and update

```bash
# Prepare metadata and every replicon for an offline prediction.
coprarna genomes NC_000913 NC_000964 NC_002516 \
  --genome-db /data/coprarna-genomes --download

# Refresh only the selected metadata when checking for newly available versions.
coprarna genomes NC_000913 --genome-db /data/coprarna-genomes --refresh

# Refresh the annotations/genomes as well, including unchanged sequence versions.
coprarna genomes NC_000913 --genome-db /data/coprarna-genomes --refresh --download
```

Use the same query IDs during prefetch and prediction. Resolutions are stored as
`metadata/<query>.json`; compressed genomes are `<accession.version>.gb.gz`.
Cached metadata stays fixed until explicitly refreshed. NCBI can change
annotations without changing the sequence accession version, so retain a separate
snapshot of a database when exact reproducibility is required. Refresh an idle
database, then use a fresh job directory for the updated prediction.

Jobs use compressed genome paths directly; they do not copy or decompress genomes
into every job. Downloads are streamed, checked for truncation and published by
atomic replacement. Per-accession locks serialize concurrent cache misses.
Failed downloads or refreshes leave the previous cache entry intact. Requests
have timeouts, bounded retries and per-process pacing. Multiple independent
processes still share NCBI's service rate limits.

Existing `.gbk.gz` and `.gbff.gz` files are also accepted. Supplied genome filenames
must identify the requested sequence version. Plain `.gb`, `.gbk` or `.gbff`
files placed in the job directory are preserved and compressed into the database
once. An offline run needs both resolution metadata and all associated genomes;
providing only a chromosome file cannot establish its assembly membership.

Without `--genome-db`, each job uses its own `genomes/` directory. Share one database
between jobs to avoid repeated downloads. `genome_metadata.json` snapshots each
job's chosen assemblies and replicons; `genome_paths.txt` records the actual files.

### Bundled reference genomes

[`core_genome_archive/`](core_genome_archive/) contains 17 refreshed RefSeq records
in directly readable gzip format. [`manifest.json`](core_genome_archive/manifest.json)
records retrieval time, exact accessions, annotation dates, lengths, CDS counts,
source URLs and SHA-256 hashes. Its metadata can seed a shared database:

```bash
coprarna genomes NC_000913 --genome-db "$PWD/core_genome_archive" --download
```

The archive contains the 17 reference replicons, not every replicon of all their
assemblies; `--download` fetches any missing associated chromosomes/plasmids.
Maintainers can refresh these references from the repository root with
`python -m tools.refresh_core_genomes`. Downloads are staged and validated before
replacing the existing reference files.

## Prediction options

### IntaRNA parameters

Predictions use **IntaRNAsTar**. Supply native IntaRNA settings in a parameter file
or as inline contents, never both. The file uses `longName=value` syntax without
`--` prefixes; see the [IntaRNA documentation](https://backofenlab.github.io/IntaRNA/).
For example, `/data/intarna.cfg`:

```ini
tAccW=80
tAccL=60
seedBP=6
```

```bash
coprarna run --workdir /data/jobs/file-settings --srnaseq /data/my_sRNAs.fa \
  --intarna-params /data/intarna.cfg --genome-db /data/coprarna-genomes
coprarna run --workdir /data/jobs/inline-settings --srnaseq /data/my_sRNAs.fa \
  --intarna-params-text $'tAccW=80\ntAccL=60\nseedBP=6\n' \
  --genome-db /data/coprarna-genomes
```

Every native prediction, including web exports, receives these settings.
`IntaRNA_parameters.txt` keeps a reproducible snapshot; temporary inline files are
removed after success or failure. CopraRNA owns the native input/output paths,
CSV format and required columns. It requests two interactions unless `outNumber`
is supplied. Native options are not duplicated as CopraRNA flags.

Genomic targets include extra sequence on both sides for accessibility calculation:
the maximum target base-pair span `tAccL` (100 nt with IntaRNAsTar defaults), or the
shared `accL` setting. With `tAccL=0`, the flank uses `tAccW` (150 nt by default,
also configurable through `accW`). If both are zero, the context covers the whole
available linear replicon, or one full turn on each side of a circular target;
this unlimited setting can be expensive for genome-wide prediction. Circular
context wraps around the origin; linear context stops at the physical chromosome
ends. The same rules apply to 5′ UTR, 3′ UTR and spliced CDS targets, on both strands.
The existing exclusion of truncated UTR targets shorter than 150 nt is unchanged.

CopraRNA sets IntaRNA's `--tRegion` to allow interactions only within the requested
region, and controls `tIdxPos0`/`qIdxPos0` so report coordinates start at 1. Targets
with different allowed intervals run in separate native batches. Accessibility is
computed on the extended sequence, including for web exports. FASTAs, interaction
CSVs, perfect-helix filtering, comparative statistics, alignments and plots use only
the requested sequence and coordinates; `E_norm` uses the requested target length.
Accessibility energies and probabilities retain the extra context's contribution,
so predicted interactions and rankings can change from earlier runs.

Each target FASTA has a `.fa.context.json` sidecar containing the extended sequence
and its 0-based, half-open requested interval. These sidecars participate in cache
signatures and are retained in `IntaRNA/` after cleanup; temporary native inputs are
removed after success or failure. Existing predictions without this context must
be recomputed. Regenerate predictions before calling
`webserver.export_prediction()` independently on an old job.

Before fitting energy distributions, CopraRNA discards gap-free, stacking-only
perfect helices spanning at least 90% of either RNA's effective length. Effective
length is the smaller of sequence length and a positive `tAccW`/`qAccW`; zero
uses the full sequence. These interactions can be artificial self hits. The same
filter applies before choosing web-report interactions. Each `.filter.json`
records removed IDs/coordinates and counts. GEV fitting requires at least three
retained optimal interactions with nonconstant energies; insufficient data stops
with a diagnostic. Cached predictions require matching input contents, parameters,
IntaRNA version and filtering policy.

### CopraRNA options

`coprarna run --help` lists accepted flags. Defaults:

| Option | Default | Meaning |
| --- | --- | --- |
| `--srnaseq` | `input_sRNA.fa` | Input relative to `--workdir`, or absolute |
| `--workdir` | Current directory | Isolated job directory |
| `--genome-db` | `<workdir>/genomes` | Shared compressed database; relative paths are job-relative |
| `--ntup`, `--ntdown` | 200, 100 | Nucleotides upstream/downstream of the selected target boundary |
| `--region` | `5utr` | `5utr` around CDS starts, `3utr` around stops, or `cds` with flanks |
| `--cores` | 1 | Native parallel work per job |
| `--nooi` | Off | Balanced ranking across organisms instead of focusing on the first |
| `--cons` | 0 | 0: no consensus weighting; 1: organism-of-interest consensus; 2: overall consensus |
| `--ooifilt` | 0 | Prefer organism-of-interest interactions with p-value at most this threshold (0 disables) |
| `--root` | 1 | Root applied to phylogenetic weights |
| `--topcount` | 200 | Top result count and extended plot count |
| `--websrv` | Off | Web tables, mapping files and SVG thumbnails |
| `--enrich` | 0 | Request DAVID enrichment for this many top predictions |
| `--noclean` | Off | Retain intermediate files for inspection |
| `--verbose` | Off | Additional progress logging |

`--cons 1` requires an organism of interest and cannot be combined with `--nooi`.
For a balanced, conservation-aware ranking:

```bash
coprarna run --workdir /data/jobs/balanced --srnaseq /data/my_sRNAs.fa \
  --nooi --cons 2 --genome-db /data/coprarna-genomes
```

## Outputs and interpretation

| File or directory | Contents |
| --- | --- |
| `CopraRNA_result.csv` | Top ranked target predictions |
| `CopraRNA_result_all.csv` | Complete comparative ranking |
| `CopraRNA_option_file.txt` | Effective job settings and software versions |
| `IntaRNA_parameters.txt` | Supplied native parameter snapshot |
| `genome_metadata.json`, `genome_paths.txt` | Assembly identities and genome paths |
| `IntaRNA/` | Gzipped native interactions and filtering audit JSON |
| `FASTA/` | Normalized input, `resolved_sRNA.fa` with primary-chromosome IDs, and target/query FASTAs |
| `Phylogeny/` | 16S alignments, distances, tree and `compatible.treefile.svg` |
| `Regions_plots/` | Target/sRNA interaction-region plots in PDF and SVG |
| `*_conservation_heatmap.pdf`, `.svg` | Comparative interaction/conservation heatmaps |
| `evo_alignments.zip` | Jalview alignments and interaction/consensus annotations |
| `weights.warning` | Organisms with unusually large phylogenetic weights |
| `README.txt` | Report archive guide |

Result columns include FDR-adjusted and comparative p-values, per-organism locus
annotations and interaction coordinates, additional homologs and retained report
fields. Coordinates follow IntaRNA's 1-based conventions relative to the extracted
sequences. Final result and web tables write unavailable values as empty cells
rather than `NA`; a missing comparative score is not evidence of a significant
interaction. Ranking identifies candidates for biological follow-up, not validated
regulation. Review the 16S tree and weight warnings when an outlier dominates.

All scientific diagrams are vector PDF/SVG. Extended PDFs retain every page, and
the matching SVG stacks those pages. `--websrv` adds `coprarna_websrv_table.csv`,
interaction details and SVG thumbnails. `--enrich N` writes optional functional
reports under `Enrichment/`. DAVID requires a reachable compatible SOAP service;
offline tests validate
request contracts, not live service availability. The optional interactive HTML
enrichment report loads upstream JavaScript when opened. Core prediction does not
need a browser or DAVID.

## Independent batches and helpers

Create separate, non-nested job directories and place `jobs.json` beside them:

```json
[
  {"id":"sample-a","helper":"predict","args":["input_sRNA.fa"],"workdir":"sample-a","options":{"cores":4,"genome_db":"/data/coprarna-genomes"}},
  {"id":"sample-b","helper":"predict","args":["input_sRNA.fa"],"workdir":"sample-b","options":{"cores":4,"cons":2,"genome_db":"/data/coprarna-genomes"}}
]
```

```bash
coprarna batch jobs.json --jobs 2 > batch-results.json
coprarna helpers
coprarna helper --workdir /data/job genome-regions NC_000913.3.gb.gz 200 100 5utr
```

`--jobs` bounds concurrent jobs; each prediction's `cores` bounds native work.
Preflight checks the complete manifest before starting. Results retain manifest
order, and other jobs continue if one fails. Helper names and positional arguments
come from `coprarna helpers`. Helpers write their documented intermediate files in
`--workdir` and return textual output on stdout. Diagnostics go to stderr.
Exit codes are 0 for success, 1 for failed batch jobs, 2 for invalid invocation or
setup, and 130 for interruption.

The optional `rasterize` helper needs `python -m pip install '.[web]'` followed by
`python -m playwright install chromium --only-shell`. It is independent of normal
scientific PDF/SVG generation.

## Upgrading from 2.1

Install into a fresh environment and start new job directories. The entry point is
`coprarna run`, replacing the original controller and script collection. Original
script-name aliases, helper defect-emulation flags, migration inventory commands
and organism-table builders have been removed. Version 3 uses only the CopraRNA2
comparative algorithm. Native accessibility/window settings belong in an IntaRNA
parameter file. The old source and documentation remain available in Git history.

The `CopraRNA_available_organisms.txt` and `kegg2refseqnew.csv` files are no longer
inputs, generated outputs or package resources. Use `coprarna genomes` to inspect,
prefetch or refresh selected assemblies instead. Preserve old databases if they
are needed to reproduce earlier results. The IntaRNAsTar defaults, helix filtering,
corrected algorithms and updated genome annotations can change rankings relative
to version 2.

## Development and validation

```bash
python -m pip install '.[test,build]'
python -m pytest -q
python -m build
COPRARNA_PIPELINE_TESTS=1 python -m pytest tests/test_workflow.py
COPRARNA_REFERENCE_TESTS=1 python -m pytest tests/test_reference_genomes.py
COPRARNA_BROWSER_TESTS=1 python -m pytest tests/test_adapters.py
```

CI checks Python 3.11–3.13, wheel/sdist builds, installed entry points and the native
scientific workflow. Tests retain independently captured statistical/helper
expectations: GEV fits, nonlinear correlation fits, rank aggregation and table
transformations. R-data matrix files are regression data read semantically by
Python; no R interpreter is invoked. Captured outputs have checksums and must not
be regenerated from the implementation being tested. Explicit contracts cover
requested behavior changes, NCBI pagination/cache failures and genome identity.
Real native tests exercise BLAST+/DomClust, IntaRNAsTar, file/inline parameter
equivalence, shared gzip genomes, input preservation and vector outputs.

Nonlinear fitting amplifies BLAS rounding differences. Validated Linux
x86-64/AVX2 runs use OpenBLAS's Haswell kernels. CopraRNA selects this before NumPy
loads when AVX2 is available, respecting explicit caller settings. If importing
NumPy first, launch with `OPENBLAS_CORETYPE=Haswell`. Exact byte parity on other
BLAS/CPU configurations is not certified.

[`AGENTS.md`](AGENTS.md) describes the module boundaries, documentation requirements
and validation policy. The project is [MIT licensed](LICENSE).

## Citation

Please cite the underlying methods when using CopraRNA:

- Wright PR, Richter AS, Papenfort K, Mann M, Vogel J, Hess WR, Backofen R, Georg J.
  [Comparative genomics boosts target prediction for bacterial small RNAs](https://doi.org/10.1073/pnas.1303248110).
  *PNAS* 110, E3487–E3496 (2013).
- Wright PR, Georg J, Mann M, Sorescu DA, Richter AS, Lott S, Kleinkauf R, Hess WR,
  Backofen R. [CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains](https://doi.org/10.1093/nar/gku359).
  *Nucleic Acids Research* 42, W119–W123 (2014).

For ad hoc use, the [Freiburg CopraRNA webserver](https://rna.informatik.uni-freiburg.de/CopraRNA/)
provides a hosted interface; its deployed version can differ from this release.
