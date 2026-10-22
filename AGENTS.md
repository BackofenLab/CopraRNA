# CopraRNA 3 architecture and release policy

[Issue #4](https://github.com/BackofenLab/CopraRNA/issues/4) prepares version 3.0.0
from `issue-1-python3-bugfixes`. Its explicit removal of the 2.1 code base supersedes
the earlier instruction to retain original source files. Historical sources remain
in Git history; production and tests must never invoke them. Keep substantial
review changes in separate commits.

## Architecture

- `pipeline.py`, `cli.py`: configuration, isolated batches, cleanup and genome CLI.
- `organisms.py`: on-demand NCBI assembly resolution, versioned metadata, cache
  validation and duplicate-assembly rejection; no static organism allow-list.
- `homology.py`, `genome.py`, `ncbi.py`, `annotation.py`: shared compressed genomes,
  streaming GenBank parsing, short job-local native species codes and clustering.
- `intarna.py`, `parameters.py`: IntaRNA >=3 / IntaRNAsTar, native parameter files,
  gzipped interactions and removal of near-full-length perfect helices.
- `energy.py`, `combination.py`, `numerics.py`, `statistics.py`: validated GEV and
  comparative statistics; preserve the validated NumPy/SciPy and BLAS kernel pins.
- `conservation.py`, `alignments.py`: consensus matrices and Jalview exports.
- `plots.py`, `phylogeny.py`: PDF/SVG scientific figures and tree rendering.
- `david.py`, `enrichment.py`, `webserver.py`: optional service and web exports.
- `io.py`: shared compressed/plain I/O and protected input/output namespaces.
- `helpers.py`: active Python helper contracts; no original-script aliases or
  defect-emulation mode. `resources.py` provides current report resources only.

## Documentation

Document each new/changed public function with its purpose, input/output contract,
file side effects and failure behavior where relevant. For scientific algorithms,
explain coordinate systems, matrix dimensions, statistical conventions and retained
oracle differences. Comment non-obvious numerical choices, native interfaces and
cache publication/lifetime rules. Record release-facing changes and rationale in
`ChangeLog.md`; keep `README.md` consistent with CLI and packaging.

## Validation

Retain exact independently captured scientific fixtures for unchanged behavior.
Their provenance and checksums live in `tests/fixtures`; do not regenerate expected
numbers from the Python implementation. Compare R-data objects semantically.
Test changed behavior against explicit contracts and real native integrations.
No legacy source/interpreter environments or generated availability catalogs belong
in the release. The current source version, package metadata and option-file
version must agree.

Native integration uses `COPRARNA_PIPELINE_TESTS=1`; reference-genome extraction
uses `COPRARNA_REFERENCE_TESTS=1`; optional browser integration uses
`COPRARNA_BROWSER_TESTS=1`. CI checks Python 3.11–3.13, builds and installed wheels.
Refresh the 17 NCBI reference records with `python -m tools.refresh_core_genomes`;
verify exact identities, checksums and annotation dates in the archive manifest.
Live external-service availability is not certified by offline request tests.
Check new dependency/UI releases before changing pins.

## Invariants

Preserve the user's original FASTA and native parameter file, including inputs
named `input_sRNA.fa`. When adding outputs, update `io._JOB_OUTPUTS`; preserve
cleanup destinations that resolve to original inputs. Every independent native
entry point retains its version check and detected-version cache signature.

Keep genomes compressed and reuse shared paths directly. Cache metadata pins exact
assembly/sequence versions until explicitly refreshed. Publish only completed
validated downloads/resolutions, under accession locks; failed refreshes preserve
prior entries. A chromosome and plasmid from the same assembly must never count
as independent comparative organisms. Scientific report labels are unversioned;
job metadata and cache paths preserve exact versions. SVG/PDF remain primary
figure formats; no PostScript/Ghostscript dependency is permitted.
