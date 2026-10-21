# CopraRNA Python 3 migration

This branch starts the behavior-preserving rewrite requested in
[BackofenLab/CopraRNA#1](https://github.com/BackofenLab/CopraRNA/issues/1).
It currently provides **22 validated preprocessing/export helpers out of 41 legacy
scripts**. It is an installable Python package with a batch interface. **It cannot
yet run a complete CopraRNA prediction.** Continue using the original
`CopraRNA2.pl` workflow for whole-genome predictions.

The original scripts, README, reference genomes and organism tables remain byte
identical to upstream commit `ee0c971`. A checksum test enforces that invariant.
The Python helpers never execute Perl, R, Bash, JavaScript or Python 2.

## Install the helper package

Python 3.11 or later and Biopython 1.84 or later are required. No external
scientific executable is needed by the currently migrated helpers.

```bash
python3 -m pip install .
coprarna --version
coprarna helpers
coprarna status
```

For development:

```bash
python3 -m pip install '.[test,build]'
python3 -m pytest
```

Alternatively, create the helper environment and install without changing its
resolved dependencies:

```bash
conda env create -f environment.yml
conda activate coprarna-python3
python -m pip install --no-deps .
```

`environment.yml` pins the helper runtime and development tools. It is not a
Bioconda recipe or an explicit platform lockfile. `environment-legacy.yml` records
historical interpreter and scientific-tool requirements for reproducing the old
workflow. That separate historical environment has not been solved/validated;
old channels/builds may require archival packages. Combining its old ViennaRNA
and MAFFT builds with a current Python environment fails dependency resolution.

Future whole-workflow packaging also needs the legacy-compatible IntaRNA,
ViennaRNA, BLAST legacy, DomClust, EMBOSS `distmat`, PHYLIP `fneighbor`, MAFFT and
Clustal Omega executables, plus replacements for the old plotting and SOAP
clients. Updating any scientific-tool version requires its own parity validation.

## Run a helper

Each helper keeps the legacy argument order and output filenames. Files are
resolved relative to `--workdir`; stdout contains only the helper's data.
`coprarna helpers` lists every available command and positional argument.
Original script basenames are accepted as aliases.

```bash
coprarna helper --workdir /path/to/job format-fasta input.fa > formatted.fa
coprarna helper --workdir /path/to/job genome-regions NC_000913.gb 200 100 5utr > targets.fa
coprarna helper --workdir /path/to/job compute-weights compatible.treefile compatible.fneighbor > zscore.weight
coprarna helper --workdir /path/to/job adjust-pvalues pvalues
coprarna helper --workdir /path/to/job prep_cop2.R opt_tags.clustered
```

The first two redirections above are ordinary shell redirections, so their output
paths are relative to your shell directory. The `adjust-pvalues` helper writes
`padj.csv` inside its working directory, as the legacy script does. Helpers that
write generated files preserve the legacy overwrite behavior. Use separate job
directories with copies of the inputs you want to process.

## Batch processing

Write `jobs.json` beside your separate job directories:

```json
[
  {"id": "sample-a", "helper": "genome-regions", "args": ["NC_000913.gb", "200", "100", "5utr"], "workdir": "sample-a", "stdout": "targets.fa"},
  {"id": "sample-b", "helper": "genome-regions", "args": ["NC_003197.gb", "200", "100", "5utr"], "workdir": "sample-b", "stdout": "targets.fa"}
]
```

```bash
coprarna batch jobs.json --jobs 4 > batch-results.json
```

Job directories must already exist and must be distinct and non-nested. Relative
job paths are resolved against the manifest's directory. The optional `stdout`
path stays inside the job directory and must not already exist. Without it,
helper output appears in that job's JSON result. Give each job its own inputs;
absolute input arguments can otherwise refer to shared files that some legacy
helpers modify.

Preflight validates the entire manifest before any job runs. Independent jobs
continue if another fails. Results retain manifest order and include the job ID,
status and either output or an error. Diagnostics go to stderr. Exit codes: `0`
for success, `1` when batch jobs fail, `2` for invalid invocation/setup, and `130`
for an interrupted command. `--jobs` controls batch concurrency; it does not
change any helper's algorithm.

## Regression tests

`tests/fixtures/*/expected` contains outputs captured from the original scripts
before their counterparts were implemented. `case.json` records interpreter,
script, arguments and generated files. Tests compare every captured stdout/file
byte. Only the bodies of explicitly marked Perl hash-order outputs are sorted;
their original raw order is nondeterministic. No numeric tolerance is used by
these golden comparisons.

Offline Python tests run without Perl or R. Legacy-oracle tests require the
original modules and automatically skip if the runtimes are unavailable:

```bash
export COPRARNA_LEGACY_PERL=/path/to/oracle/bin/perl
export COPRARNA_LEGACY_RSCRIPT=/path/to/oracle/bin/Rscript
python3 -m pytest
```

The checked-in outputs were captured using Perl 5.26.2 with BioPerl and R 4.4.3.
The local oracle versions differ from the historical README's R 3.2.2, so those
historical numerical versions have not been certified. There are no upstream
example sRNA inputs or full-run expected outputs in this repository.

The optional real-genome check decompresses the 17 bundled genomes into temporary
directories and compares all three extraction modes directly against Perl:

```bash
COPRARNA_REFERENCE_TESTS=1 python3 -m pytest tests/test_reference_genomes.py
```

To intentionally rebuild small synthetic inputs and then refresh their legacy
outputs after reviewing changes:

```bash
python3 tools/build_test_fixtures.py
python3 tools/capture_legacy.py
```

The GEV fixture is an **oracle baseline only**: its Python port does not exist
yet. Passing helper tests does not establish end-to-end prediction parity.

## Script mapping and remaining work

Each original script has a neighboring `<script-name>.md` mapping document.
`coprarna/migration.json` records implemented targets, commands, fixtures and
planned targets for unimplemented scripts. `coprarna status` prints this inventory.
`AGENTS.md` is the architectural ledger; regenerate mapping documents with
`python3 tools/update_migration_inventory.py` after adding validated helpers.

The remaining stages are the full pipeline and NCBI retrieval, IntaRNA and
alignment orchestration, generalized-extreme-value fitting, correlated CopraRNA1
and CopraRNA2 joining, consensus/conservation calculations, R graphics and Rdata
compatibility, DAVID SOAP clients, webserver rendering and organism-list updates.
These stages must be translated and validated before the rewrite PR is ready.
In particular, a generic SciPy fit must not silently replace R's `optim`/`nls`,
and NumPy sampling must not be assumed to reproduce R's random sequence.

Known legacy issues are recorded in `AGENTS.md`. The incorrect column-35 rename
in the cluster-size filter is deliberately preserved, as are the region parser's
first-record sequence behavior and its four-base complement handling. Corrections
belong in a separate branch and follow-up PR after the parity rewrite. No legacy
bug fixes are included here.

The original [README.md](README.md) contains biological context, citations and
instructions for the unchanged legacy workflow. The project remains MIT licensed.
