# CopraRNA with Python 3

The complete parallel rewrite requested in [issue 1](https://github.com/BackofenLab/CopraRNA/issues/1)
provides Python counterparts for all 41 legacy scripts. It includes comparative
prediction, CopraRNA1 and CopraRNA2 statistics, consensus mapping, plots, enrichment,
webserver tables, organism-list maintenance and independent prediction batches.
The original sources, README, reference genomes and organism tables remain
unchanged; a checksum test verifies every original tracked file.

## Installation

Python 3.11 or newer is required. Create the full scientific environment, then
install the package:

```bash
conda env create -f environment.yml
conda activate coprarna-python3
python -m pip install --no-deps .
coprarna --version
```

`environment.yml` pins the validated Linux runtime: Biopython, NumPy, SciPy,
pandas, Matplotlib, rdata, suds-community, IntaRNA 3.4.1, ViennaRNA 2.7.2,
BLAST legacy 2.2.26, DomClust 1.2.8, EMBOSS 6.6.0, PHYLIP fneighbor,
MAFFT 7.526, Clustal Omega 1.2.4 and Ghostscript. This specification was installed
in a fresh environment. It prepares dependency metadata for a future Bioconda
recipe; it is not a cross-platform lockfile or a published Bioconda package.

For Python-only helper/statistical use, `python -m pip install .` is sufficient.
Install native scientific executables separately for complete predictions.
`COPRARNA_TOOL_PATH` can prepend their directory to subprocess PATH. Production
controllers do not execute the legacy Perl, R, Bash, Python or JavaScript scripts.
Native tools can have their own transitive interpreter dependencies.

For enrichment graphics, install Chromium once:

```bash
python -m pip install '.[web]'
python -m playwright install chromium --only-shell
```

DAVID enrichment requires a reachable compatible SOAP service and its legacy
registered account. Requests preserve the original endpoint, category settings,
background population and account. NCBI retrieval needs network access unless
all required GenBank replicons already exist in the working directory. These
services are tested with offline request doubles; current live service availability
is not promised. The unchanged HTML reports retain their upstream web dependencies.

## Predict targets

Prepare at least three homologous sRNAs in FASTA. Each header must be an organism
RefSeq accession present in the bundled availability table. Put this file in a
separate job directory. Cached `<accession>.gb` files avoid NCBI downloads.

```bash
coprarna run --workdir /data/job-a --srnaseq input_sRNA.fa --cores 8
coprarna run --workdir /data/job-b --srnaseq input_sRNA.fa --nooi --cons 2
coprarna run --workdir /data/job-c --srnaseq input_sRNA.fa --cop1 --seed 731
coprarna run --help
```

Options retain the legacy defaults and accept both `--name` and `-name` spellings:
`ntup`, `ntdown`, `region` (5utr/3utr/cds), `rcsize`, `cores`, `winsize`,
`maxbpdist`, `cop1`, `verbose`, `websrv`, `nooi`, `topcount`, `root`, `enrich`,
`noclean`, `cons` and `ooifilt`. `--seed` makes CopraRNA1's missing-value sampling
reproducible. An unspecified seed retains stochastic behavior.

The workflow retrieves/caches genomes, extracts proteins and target regions,
builds homology clusters, computes IntaRNA interactions, estimates GEV p-values,
weights organisms using 16S phylogeny, combines predictions and maps conservation.
Outputs include `CopraRNA_result.csv`, `CopraRNA_result_all.csv`, evolutionary
heatmaps, interaction-region plots and a Jalview alignment archive. `--noclean`
retains intermediate files. `--websrv` adds the webserver interaction table;
`--enrich N` requests functional enrichment of the first N predictions.

This rewrite branch preserves documented legacy defects, including DomClust's
obsolete table-output mode and failures when plotting lowercase cds or fewer
than the requested predictions. An existing `cluster.tab` uses the original
supported imported-cluster path. Corrections are delivered separately on the
stacked bug-fix branch. See `AGENTS.md` for the exact compatibility decisions.
Use fresh job directories: legacy append/overwrite behavior is preserved.

## Independent batches

Place `jobs.json` beside separate, non-nested job directories:

```json
[
  {"id":"sample-a","helper":"predict","args":["input_sRNA.fa"],"workdir":"sample-a","options":{"cores":4,"topcount":100}},
  {"id":"sample-b","helper":"predict","args":["input_sRNA.fa"],"workdir":"sample-b","options":{"cores":4,"cons":2}}
]
```

```bash
coprarna batch jobs.json --jobs 2 > batch-results.json
```

`--jobs` limits concurrent jobs; each prediction's `cores` limits its native
parallel work. Preflight validates the complete manifest before starting any job.
Results retain manifest order, and independent jobs continue after another fails.
Diagnostics go to stderr. Exit codes are 0 for success, 1 for failed batch jobs,
2 for invalid invocation/setup, and 130 for interruption.

Helpers can also be batched by replacing `predict` with a helper name, supplying
its positional arguments and omitting `options`. Optional `stdout` writes a
new file inside that job directory. `coprarna helpers` lists arguments; original
script basenames are accepted as aliases.

```bash
coprarna helper --workdir /data/job-a genome-regions NC_000913.gb 200 100 5utr
coprarna helper --workdir /data/job-a fit-energy
coprarna status
```

## Validation and development

```bash
python -m pip install '.[test,build]'
python -m pytest -q
python -m build
```

The committed small fixtures contain independently captured legacy stdout and
files. Python comparisons use exact bytes, except explicitly unordered Perl hash
bodies. The 120-cluster statistical fixtures cover missing homologs, R's seeded
random stream, GEV fitting, nonlinear correlation fits and rank aggregation.

`tests/fixtures/workflow/legacy.tar.gz` captures an offline three-genome,
30-cluster synthetic run using the imported-cluster path. Tests compare all
prediction CSVs, target alignments, 203 Jalview files, plot coordinates, heatmap
scores/labels/order and persisted matrices. R-data comparisons check object
contents, dimensions and names because serialization compression/header bytes
can differ. Matplotlib and Chromium replace R graphics and PhantomJS; rendered
PDF/PS/PNG bytes and typography are not identical across renderers. PDF and
PostScript support multiple pages. External-tool versions are recorded in fixture
provenance; historical R 3.2.2 and IntaRNA 2.2.0 numerical equivalence is not certified.

```bash
COPRARNA_PIPELINE_TESTS=1 python -m pytest tests/test_workflow.py
COPRARNA_BROWSER_TESTS=1 python -m pytest tests/test_adapters.py
COPRARNA_REFERENCE_TESTS=1 python -m pytest tests/test_reference_genomes.py
```

The last command additionally requires the old Perl dependencies and compares
all three regions for each of the 17 bundled real genomes. To rerun original
oracles, set `COPRARNA_LEGACY_PERL` and `COPRARNA_LEGACY_RSCRIPT` to suitable
interpreters. Small fixtures were captured with Perl 5.26.2/R 4.4.3; full-workflow
and report fixtures used Perl 5.32.1/R 4.5.2. `environment-legacy.yml` records the
historical requirements separately and is not the production environment.

Each legacy script has an adjacent `.md` mapping. `coprarna/migration.json`
identifies its Python function, command and regression coverage. Regenerate these
with `python tools/update_migration_inventory.py`. `AGENTS.md` records architecture,
validation and deferred defects. The original [README.md](README.md) retains
biological background and citations. The project remains MIT licensed.

The nonlinear fit's finite differences amplify BLAS rounding differences. The
validated Linux x86-64/AVX2 runtime pins NumPy 2.4.3, SciPy 1.17.1 and OpenBLAS's
Haswell reduction kernels (also used on Zen). CopraRNA selects this kernel before
NumPy loads when AVX2 is available. If importing NumPy first, launch with
`OPENBLAS_CORETYPE=Haswell`; an explicit caller setting is respected. Exact byte
parity on other BLAS/CPU configurations is not certified.
