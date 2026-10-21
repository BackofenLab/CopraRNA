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
