# `homology_intaRNA.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.homology.run_homology`.

Command: `coprarna helper homology` followed by `srna up down region accession`.

Regression fixtures: `tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
