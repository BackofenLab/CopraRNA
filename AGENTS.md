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
