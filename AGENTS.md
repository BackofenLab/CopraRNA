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
