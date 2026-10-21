# CopraRNA Python 3 migration ledger

## Scope and invariants

Issue: https://github.com/BackofenLab/CopraRNA/issues/1.

Create an installable Python 3 implementation alongside the unchanged legacy
sources. Preserve algorithms, heuristics, numerical behavior, ordering and output
formats. Record defects here; do not silently fix them in the parity rewrite.
Keep each major migration step in an atomic commit. Bug fixes belong on a separate
branch and pull request after the parity implementation.

## Architecture and module mapping

- `CopraRNA2.pl`: pipeline and legacy CLI; destination `coprarna/pipeline.py` and
  `coprarna/cli.py`.
- `coprarna_aux/*.pl`: parsing, clustering, annotation and formatting helpers;
  destinations under `coprarna/` with per-script mapping documents alongside the
  originals.
- `coprarna_aux/*.R`: statistical combination, filtering and plotting; Python
  modules must be validated against legacy R outputs before claiming parity.
- `coprarna_aux/*.py`: DAVID service clients; Python 3 service adapters.
- `update_kegg2refseq/*`: organism mapping maintenance; Python maintenance module.
- `tests/`: offline fixtures, exact-output legacy regression comparisons and
  optional integration checks requiring external scientific tools.
- `environment.yml`, `pyproject.toml`, `README.2.md`: reproducible dependencies,
  package metadata and documentation for the new interface.

External scientific executables remain external dependencies; the new runtime
must not invoke Perl, R, Bash or Python 2 scripts to implement migrated logic.

## Workflow

Inspect each legacy component before porting. Capture deterministic legacy outputs
and compare bytes where applicable. Document missing runtimes, network resources,
fixtures and numerical validation explicitly. Never claim full workflow parity
from helper-level checks. Preserve all tracked original source bytes.

## Progress

- Repository cloned into `/home/alex/CodexProjects/CopraRNA`.
- Phase 0: architectural ledger initialized before implementation.
- Phases 1–6: pending repository audit and baseline validation.

## Legacy defects and compatibility decisions

Pending audit. Fixes require separate tests and a separate branch/PR.

### Phase 1 audit

No upstream tests, sRNA examples or expected result directories are tracked.
The 17 compressed GenBank reference genomes can support real-data parser tests.
`tests/support/legacy.py` captures stdout and declared files in isolated temporary
working directories; golden fixtures are versioned. Oracle runtimes may be supplied
through `COPRARNA_LEGACY_PERL` and `COPRARNA_LEGACY_RSCRIPT`.
The local reusable oracle runtimes are Perl 5.26.2 (virsorter1 environment) and
R 4.4.3 (virfinder environment). These differ from the historical versions in the
README; numerical parity claims must specify the oracle actually used.

The migration environment is specified separately from the legacy comparison
environment. Conda specifications are compatibility targets, not solved locks.

Known parity hazards discovered during audit:
- Perl hash iteration randomizes output row order in several helpers. Fixtures
  mark those files explicitly; comparisons sort only the unordered body rows.
- `join_pvals_coprarna1.R` samples missing values without a fixed random seed.
- `remove_clusters_under_percantage.R` renames column 35 rather than 36.
- `parse_region_from_genome.pl` always extracts sequence from the first GenBank
  record, including when iterating features of later records.

Phase 1 validation: 28 captured legacy fixture cases plus the original-source
checksum test pass (29 tests). All golden files were captured before porting.
The fixtures include a full GEV fit oracle; this does not establish Python GEV
parity. Oracle commands used the explicit runtime paths described above.

### Phase 2: sequence, table and phylogeny helpers

Implemented Python ports in `fasta.py`, `genome.py`, `tables.py`, `annotation.py`,
`phylogeny.py` and `statistics.py`. The 28 migrated synthetic fixture cases match
legacy outputs (byte comparisons, except explicitly unordered Perl hash bodies).

Real-genome validation exposed a translation error in compound locations spanning
the circular origin. BioPerl keeps the location-list order; Biopython reverses a
complemented location list. The port now reverses only complemented lists and
preserves origin-spanning joins. This corrects the translation, not legacy logic.
The affected NC_000911 genome now matches all three regions; the other 16 genomes
already matched all three regions. Total: all 51 real-genome/region comparisons
have passed with the final location handling (16 genomes before the correction;
one affected genome afterward).

Not yet ported: GEV fitting, correlated p-value joining, consensus/conservation,
plotting, full pipeline orchestration, webserver rendering and maintenance tools.
No complete end-to-end rewrite or whole-workflow equivalence is claimed.

### Export helpers, CLI and packaging

Added byte-validated enrichment JSON and specific-target exports, the unchanged
archive README resource, and compressed copies of the unchanged organism tables.
`helpers.py` is the production dispatch shared by the CLI and parity tests.
`cli.py` provides single-helper execution and independent bounded-concurrency
batch jobs; `pyproject.toml` packages only the Python implementation and its data.
`migration.json` and per-script `.md` files distinguish actual ports from intended
future targets. The CLI does not expose an unimplemented prediction command.

The first complete environment solve failed: historical ViennaRNA 2.4.14 builds
require Python 2, and MAFFT 7.310 conflicts with strict modern channel priority.
`environment.yml` therefore specifies only the actual helper dependencies;
whole-workflow historical tools are recorded separately in the unvalidated
`environment-legacy.yml`. A future complete Bioconda recipe remains outstanding.

Additional legacy defects retained/documented:
- Heatmap JSON is assembled without escaping feature/source strings.
- `get_CDS_from_gbk.pl` searches the mapping using the input filename minus three
  characters, so passing a path can prevent organism-code lookup.
- The region parser complements only A/C/G/T in flanking regions, not all IUPAC.
- `disentangle_sorted_intarna_CSV.pl` keeps the last suboptimal row if more than
  two interactions exist for a target.

Full rewrite completion and the separate bug-fix PR are still pending. Do not use
closing issue keywords in a draft PR for this incremental implementation.

Final validation progress:
- Full final real-genome rerun: 17/17 genomes, all three region modes passed.
- Small-fixture/oracle/CLI suite: 78 tests passed; 17 real-genome tests are
  intentionally opt-in in the default run (verified separately above).
- Additional R-derived numeric fixtures verify scientific notation and numeric
  type conversion, preventing text-format drift in the preprocessing ports.
- The revised helper Conda specification solved successfully; a fresh pinned
  environment and installed-wheel smoke checks are being validated.

Pinned-environment validation found a BioPython-version compatibility difference:
1.84 retains internal tabs in parsed FASTA sequences, while 1.86 removes them.
`fasta.format_fasta` now removes sequence whitespace explicitly, as BioPerl does.
This is a port correction, not a legacy behavioral change. All 45 Python
parity/CLI/inventory tests pass with the pinned Python 3.11 / BioPython 1.84.

Distribution checks completed:
- Fresh `environment.yml` installation succeeded on linux-64.
- The installed wheel was exercised outside the repository: all 33 migrated
  fixture cases passed through the console entry point with exact-output checks.
- The wheel contains no legacy scripts; the source archive retains every
  original tracked file with its original checksum.
- `tests/fixtures/provenance.json` records oracle and Python dependency versions.
- Ruff and `git diff --check` pass. CI covers Python 3.11–3.13, source builds and
  an installed-wheel smoke test. Historical whole-workflow environment, full
  numerical parity, the prediction CLI and the bug-fix PR remain pending.

Final complete validation: **97 tests passed** in the fresh pinned Python 3.11 /
Biopython 1.84 environment with both legacy oracle runtimes and
`COPRARNA_REFERENCE_TESTS=1` enabled. This includes all 17 bundled real genomes,
each checked in 5utr, 3utr and cds modes (51 exact-output comparisons).

CI found that Git omits the archive-README fixture's empty input directory. Added
an explicit `.gitkeep` and updated the fixture builder so clean checkouts retain
that input directory. This changes only test fixture transport, not runtime logic.

## Continuation: completing issue 1

The GEV simplex fit, CopraRNA1 seeded missing-value sampling and CopraRNA2 subset
joining are now implemented in `energy.py`, `numerics.py` and `combination.py`.
`rcompat.py` preserves table formatting and writes native R-data archives through
a pure Python serializer. New independently captured 120-cluster/four-organism
fixtures exercise full and missing homolog sets. All 73 statistical/helper oracle
and Python parity checks pass, including byte-exact CopraRNA1 and all CopraRNA2
result files. Matching the nonlinear fits required preserving R's reciprocal
scaling in LINPACK QR and the empirical-CDF multiplication order.

New legacy defect: `join_pvals_coprarna1.R` reads `h$intensities`, but histogram
objects expose `density`. The absent member yields NULL and an error score of
zero for every candidate rho. The rewrite preserves this behavior and still
consumes the same random draws; a fix belongs in the separate bug-fix branch.

The full native orchestration and reporting modules are under integration test.
A synthetic offline three-genome fixture uses the supported imported-cluster
path. The installed DomClust 1.2.8 ignores trailing options after its positional
files, so the original invocation produces list output instead of the requested
table. This invocation defect is reserved for the bug-fix branch.

### Complete workflow and reporting validation

All 41 original scripts now have Python counterparts. The full synthetic
three-genome/30-cluster workflow matches legacy prediction CSVs, weights,
phylogeny distances and target alignments byte for byte. All 203 Jalview exports
match exactly, including R's retained source row numbers and padded matrix
columns. Conservation, interaction and sparse named consensus R objects are
checked semantically; serialization headers and compression are not compared.
The seeded full CopraRNA1 workflow also matches final tables and sampled outputs.

Plot geometry, heatmap scores, annotations and average-linkage leaf order match
independently captured R data. Matplotlib replaces R rendering, so typography
and PDF/PS bytes differ. Ghostscript preserves all pages in PostScript exports.
Chromium replaces PhantomJS and is tested by real PNG/PDF rendering. SOAP clients
are compared directly with the unchanged Python source using an offline client;
NCBI request contracts and cached-file behavior are tested without live services.
The seeded organism-list fixture was captured from untouched Perl/Bash sources.

Modern IntaRNA removed old UI modes. The Python webserver adapter reproduces the
IntaRNA v2.2.0 `OutputHandlerIntaRNA1.cpp` ASCII layout from native CSV/base-pair
data, using the original Perl wrapper's scientific options. No legacy wrapper is
executed. The original CopraRNA formatter is still independently byte-validated.

The fresh full `environment.yml` install succeeded after clearing disposable task
caches. Runtime pins now include the native scientific executables and Python
statistics/reporting dependencies. Full pipeline tests are opt-in because they
require these executables; the complete offline fixture is committed for CI.
Remaining work in this task is verification/publishing and the separate bug-fix
branch required by phase 6, not additional unported components.
