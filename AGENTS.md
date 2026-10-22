# CopraRNA2 architecture and review ledger

Issue: https://github.com/BackofenLab/CopraRNA/issues/1.
Rewrite PR: https://github.com/BackofenLab/CopraRNA/pull/2.
Follow-up PR: https://github.com/BackofenLab/CopraRNA/pull/3.

## Current scope

Martin's review (https://github.com/BackofenLab/CopraRNA/pull/3#issuecomment-6324702688)
changes the follow-up to a CopraRNA2-only runtime. Its requested behavior changes
supersede historical whole-workflow parity. Remove retired-mode code, tests and
user documentation from the new implementation. Keep every original tracked
legacy source byte unchanged for historical comparison; never invoke it from
production. Keep substantial review changes in separate commits.

## Architecture

- `pipeline.py`, `cli.py`: prediction configuration, isolated batches and cleanup.
- `homology.py`, `genome.py`, `ncbi.py`, `annotation.py`: genome retrieval, shared
  compressed genome database, streaming GenBank parsing and homology clustering.
- `intarna.py`: IntaRNA >=3 / IntaRNAsTar, file or inline parameter configuration,
  gzipped interactions and removal of near-full-length perfect helices.
- `energy.py`, `combination.py`, `numerics.py`: validated GEV and comparative
  statistical algorithms. Pin the validated NumPy/SciPy and BLAS reduction kernel.
- `conservation.py`, `alignments.py`: consensus matrices and Jalview exports.
- `plots.py`, `phylogeny.py`: PDF/SVG scientific figures and tree rendering.
- `david.py`, `enrichment.py`, `webserver.py`: optional service and web exports.
- `io.py`: streaming compressed/plain file access shared by producers/consumers.
- Adjacent `.md` mapping documents and `migration.json` distinguish active and
  intentionally retired scripts; no retired helper is advertised as implemented.

## Code documentation

Document each new or changed public function with its purpose, input/output
contract, file side effects and failure behavior where relevant. For scientific
algorithms, explain coordinate systems, matrix dimensions, statistical conventions
and differences deliberately retained for oracle comparisons. Add inline comments
for non-obvious numerical choices, native-tool interface constraints and cache
publication/lifetime rules. Keep documentation accurate when changing behavior;
do not add comments that merely restate individual statements. Record release-facing
changes and their rationale in `ChangeLog.md`.

## Validation policy

Keep exact captured CopraRNA2 helper/statistical oracles where behavior is retained.
Test requested changes against explicit contracts and real native integrations.
Do not regenerate old numerical expectations from the translated implementation.
R-data objects are compared semantically; original sources have checksum tests.
SVG/PDF are primary output formats; no PostScript/Ghostscript dependency remains.
Preserve the user's original input FASTA, including when its name is input_sRNA.fa.
Keep downloaded genomes compressed and reuse shared database paths directly.

## Progress

- Previous rewrite and fixes were delivered with CI green on Python 3.11–3.13,
  exact statistical/report fixtures, 51 region comparisons on 17 reference genomes,
  native workflows, browser rendering and installed-wheel checks.
- Current review: current-only runtime, native parameter files, gzip/cache support,
  helix filtering, vector reports, documentation and current dependencies implemented.
  Fresh environment, retained numerical fixtures and native workflows pass; final
  release checks and CI are recorded in the follow-up PR.
- Live external service availability is not certified. Offline adapters retain
  request tests. New dependencies/UI versions must be checked before pin updates.


## Completeness audit — 2026-10-22

The original issue and latest review were checked together; the current review
supersedes historical defaults only where explicitly requested. Original sources
remain subject to the checksum test. No additional inline review requests exist.

| Requested behavior | Implementation and verification |
| --- | --- |
| Current comparative runtime only | `Config`, helper registry and `migration.json`; retired mode/fixtures removed, 39 active and 2 retired mappings verified by `test_inventory.py` and `test_review.py`. |
| IntaRNA >=3 and sTar default | `parameters.check_version`, native homology/prediction/web entry points; old releases rejected before side effects, plus real native workflow tests. |
| File/inline native parameters | `parameters.parameter_file`/`arguments`; mutual exclusion, exact snapshots, temporary-file lifetime and native file/inline equivalence tests. |
| Native CSV and no duplicate CLI settings | `intarna.py`, `webserver.py`, `cli.py`; gzip CSV parsing and retired CLI flag checks. |
| Current dependency interfaces | Fresh `environment.yml`, native BLAST+/DomClust integrations, exact retained numerical fixtures, Python 3.11–3.13 CI and separate optional oracle environments. |
| Compressed intermediates/genomes | `io.py`, `ncbi.py`, native `.csv.gz` producers; streamed HTTP gzip/plain handling, atomic truncation rejection, direct cache reuse and concurrent-access tests. |
| Python tree SVG | `phylogeny.render_tree`; SVG vector-path contract and native output checks. |
| PDF/SVG diagrams, no PostScript | Python plotting and enrichment; matching vector files, all extended pages, no embedded raster images and browser adapter tests. |
| Original FASTA preservation | Working copy, output-name/symlink collision preflight and guarded cleanup sources/destinations; failure-path and complete native tests. |
| Shared compressed genome database | `--genome-db`, manifest and per-accession locking; two real jobs share one cache without copying genomes. |
| Remove near-full-length perfect helices | `filter_perfect_helices` before GEV and web selection; boundary/loop cases, native helix and discarded-best-hit report regressions. |
| Function/inline documentation | Public module API docstrings, numerical/interface comments and directives above. |
| Rewrite and fix ChangeLog | `ChangeLog.md` identifies source version, architecture and each fix/change; included in sdist. |
| Original issue CLI/package/mappings | Separate rewrite and fix PRs, isolated batches, per-script sidecars, wheel/sdist and installed native prediction checks. |

When adding outputs, update `io._JOB_OUTPUTS` and preserve both original inputs
and any cleanup destination that resolves to them. Keep every independent native
entry point's version check and detected-version cache signature in place.
