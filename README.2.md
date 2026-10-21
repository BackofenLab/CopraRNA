# CopraRNA2 with Python 3

The Python implementation for [issue 1](https://github.com/BackofenLab/CopraRNA/issues/1)
provides comparative target prediction, consensus mapping, scientific plots,
functional enrichment, webserver exports and independent prediction batches.
It implements 39 active legacy-script counterparts. Two obsolete scripts are
retired in the migration inventory. Original tracked sources, reference genomes
and organism tables remain unchanged and are verified by checksums.

## Installation

Python 3.11 or newer is required. The complete Linux scientific environment uses
Python 3.12:

```bash
conda env create -f environment.yml
conda activate coprarna-python3
python -m pip install --no-deps .
coprarna --version
```

The reviewed dependency versions are listed below. IntaRNA must be version 3 or
newer. Production controllers run Python and native scientific tools; they do
not invoke the original scripts. Native packages can have transitive interpreter
dependencies. The environment specification is not a cross-platform lockfile or
a published Bioconda package.

| Dependency | Reviewed version / interface |
| --- | --- |
| [IntaRNA](https://github.com/BackofenLab/IntaRNA/releases) / ViennaRNA | 3.4.1 / 2.7.2; `IntaRNAsTar`, native gzip CSV |
| [BLAST+](https://bioconda.github.io/recipes/blast/README.html) | 2.17.0; `makeblastdb` and `blastp`, tabular format 6 |
| DomClust | 1.2.8; supported membership output mode 9 |
| EMBOSS / embassy-phylip | 6.6.0 / 3.69.660; `distmat` and `fneighbor` |
| MAFFT / Clustal Omega | 7.526 / 1.2.4 |
| [Biopython](https://pypi.org/project/biopython/) | 1.88; package also supports 1.84 on Python 3.11 |
| [NumPy](https://pypi.org/project/numpy/) / [SciPy](https://pypi.org/project/scipy/) | 2.5.4 / 1.18.1 on Python >=3.12; 2.4.3 / 1.17.1 on Python 3.11 |
| [pandas](https://pypi.org/project/pandas/) / [Matplotlib](https://pypi.org/project/matplotlib/) | 3.0.6 / 3.11.2 |
| rdata / suds-community | 1.1.0 / 1.2.0 |

Matplotlib is installed from its PyPI wheel to avoid an ICU conflict between the
current conda rendering dependencies and IntaRNA. PDF/SVG reporting uses Python;
no browser or external format converter is needed for prediction reports.
For Python-only helpers, `python -m pip install .` is sufficient. Install native
executables separately for complete predictions. `COPRARNA_TOOL_PATH` can prepend
one or more tool directories to subprocess PATH.

DAVID enrichment requires a reachable compatible SOAP service and the registered
account used by the original workflow. Offline tests check request contents;
live service availability is not certified. The optional interactive HTML report
loads its upstream JavaScript dependencies when opened in a browser.

## Predict targets

Prepare at least three homologous sRNAs in FASTA. Each record ID must be an
organism RefSeq accession in the bundled availability table. Use a separate job
directory for each run. The original input FASTA stays at its original path,
byte unchanged; a normalized working copy is stored under `FASTA/`.

```bash
coprarna run --workdir /data/job-a --srnaseq input_sRNA.fa --cores 8 \
  --genome-db /data/coprarna-genomes
coprarna run --workdir /data/job-b --srnaseq input_sRNA.fa --nooi --cons 2 \
  --genome-db /data/coprarna-genomes
coprarna run --help
```

The shared genome directory stores `<accession>.gb.gz` files. Existing `.gbk.gz`
and `.gbff.gz` files are accepted too. Jobs read these compressed files directly
without making genome copies or decompressed files in their working directories.
Downloads stream into gzip files and are published atomically after completion;
per-accession locks prevent concurrent jobs from downloading the same genome.
An existing compressed genome in the job directory is also accepted. A supplied
plain `.gb` genome is preserved and streamed into the cache once. NCBI network
access is needed only for missing replicons. `genome_paths.txt` records the paths
used by a job. Without `--genome-db`, the cache defaults to `genomes/` in the job.

### IntaRNA parameters

Predictions use the **IntaRNAsTar** personality. Supply native IntaRNA settings
through either a parameter file or its inline contents. These alternatives are
mutually exclusive. A file uses IntaRNA's `longName=value` syntax, without `--`
prefixes; see the [official parameter-file documentation](https://backofenlab.github.io/IntaRNA/).
For example, `intarna.cfg` might contain:

```ini
tAccW=80
tAccL=60
seedBP=6
```

```bash
coprarna run --workdir /data/job-a --intarna-params intarna.cfg
coprarna run --workdir /data/job-b \
  --intarna-params-text $'tAccW=80\ntAccL=60\nseedBP=6\n'
```

File paths are relative to the job directory unless absolute. Every prediction
invocation, including the web export, receives the parameter file. Inline text is
materialized temporarily and removed even if prediction fails. A persistent
`IntaRNA_parameters.txt` snapshot records the settings. CopraRNA sets its input,
output format and required CSV columns; the web export requests one interaction
per target. The primary prediction requests two interactions unless the parameter
file supplies `outNumber`. Native prediction settings are not duplicated as
CopraRNA command-line options.

Before fitting energy distributions, CopraRNA removes stacking-only, gap-free
helices spanning at least 90% of either RNA's effective length. Effective length
is the smaller of the sequence length and a positive accessibility window
(`tAccW` or `qAccW`); window zero uses the full sequence length. These interactions
can be artificial self hits. Each organism's `.filter.json` records the removed
IDs and coordinates, counts and a cache signature. The GEV fit requires at least
three retained optimal interactions and nonconstant energies; otherwise the run
stops with a diagnostic. Cached predictions are reused only when the input,
parameters, IntaRNA version and filter policy match.

Other options include `ntup`, `ntdown`, `region` (5utr/3utr/cds), `cores`,
`verbose`, `websrv`, `nooi`, `topcount`, `root`, `enrich`, `noclean`, `cons` and
`ooifilt`. Their defaults and accepted spellings are shown by `run --help`.

### Outputs

`CopraRNA_result.csv` contains the top predictions; `CopraRNA_result_all.csv`
contains all comparative results. Gzipped IntaRNA intermediates are written
natively and retained under `IntaRNA/` after cleanup. Genomes remain compressed
in the selected cache. `--noclean` keeps the other intermediate files.

`Phylogeny/compatible.treefile.svg` renders the 16S tree with Biopython and
Matplotlib. Region plots, conservation heatmaps and functional enrichment plots
are exported as **PDF and SVG**. Extended lists retain every page in the PDF and
stack all pages in one SVG canvas. Heatmaps use vector cells. `--websrv` adds the
webserver table and SVG thumbnails. `evo_alignments.zip` contains Jalview exports.
`--enrich N` requests DAVID functional enrichment of the first N predictions.

## Independent batches and helpers

Place `jobs.json` beside separate, non-nested job directories:

```json
[
  {"id":"sample-a","helper":"predict","args":["input_sRNA.fa"],"workdir":"sample-a","options":{"cores":4,"genome_db":"/data/coprarna-genomes"}},
  {"id":"sample-b","helper":"predict","args":["input_sRNA.fa"],"workdir":"sample-b","options":{"cores":4,"cons":2,"intarna_params_text":"tAccW=80\ntAccL=60\n"}}
]
```

```bash
coprarna batch jobs.json --jobs 2 > batch-results.json
coprarna helper --workdir /data/job-a genome-regions NC_000913.gb.gz 200 100 5utr
coprarna helpers
coprarna status
```

`--jobs` limits concurrent jobs; each prediction's `cores` limits native parallel
work. Preflight validates the complete manifest before starting jobs. Results
retain manifest order, and independent jobs continue after another fails.
Diagnostics go to stderr. Exit codes are 0 for success, 1 for failed batch jobs,
2 for invalid invocation/setup, and 130 for interruption.

Helpers accept original script basenames as aliases. Helper-only `--legacy`
reproduces captured defects for retained historical algorithm comparisons;
complete prediction always uses the corrected implementation. The optional
`rasterize` helper can export browser content after installing `.[web]` and
Chromium with `python -m playwright install chromium --only-shell`.

## Validation and development

```bash
python -m pip install '.[test,build]'
python -m pytest -q
python -m build
COPRARNA_PIPELINE_TESTS=1 python -m pytest tests/test_workflow.py
COPRARNA_BROWSER_TESTS=1 python -m pytest tests/test_adapters.py
COPRARNA_REFERENCE_TESTS=1 python -m pytest tests/test_reference_genomes.py
```

Retained helper/statistical fixtures use independently captured legacy output.
Byte comparisons cover GEV fitting, nonlinear correlation fits, rank aggregation
and table transformations. Workflow fixtures check consensus matrices, 203
Jalview exports, alignment inputs and plot geometry. R-data files are compared
semantically. New native integrations check the sTar workflow, BLAST+, file/inline
parameter equivalence, compressed shared genomes, input preservation and vector
outputs. Default whole-workflow predictions intentionally change with the reviewed
personality and filtering; they are not asserted equal to historical predictions.

To rerun original helper oracles, set `COPRARNA_LEGACY_PERL` and
`COPRARNA_LEGACY_RSCRIPT` to suitable interpreters. `environment-legacy.yml`
records optional interpreter requirements and the archived RobustRankAggreg
installation command separately. Its source package works with R 4.4.3 without
the obsolete ABI pins of the old conda binary. The reference-genome test
compares all three regions on all 17 bundled genomes and needs the Perl modules.
The updated archive README has its own expectation because its output changed.

Nonlinear fitting amplifies BLAS rounding differences. Validated Linux
x86-64/AVX2 runs use OpenBLAS's Haswell reduction kernels. CopraRNA selects this
kernel before NumPy loads when AVX2 is available. If importing NumPy first,
launch with `OPENBLAS_CORETYPE=Haswell`; an explicit caller setting is respected.
Exact byte parity on other BLAS/CPU configurations is not certified.

[ChangeLog.md](ChangeLog.md) records the rewrite architecture and each release-facing
fix and behavior change. Each original script has an adjacent `.md` mapping. `coprarna/migration.json`
records active functions, commands, tests and retired entries. Regenerate it with
`python tools/update_migration_inventory.py`. `AGENTS.md` records architecture and
validation policy. The original [README.md](README.md) retains biological
background and citations. The project remains MIT licensed.
