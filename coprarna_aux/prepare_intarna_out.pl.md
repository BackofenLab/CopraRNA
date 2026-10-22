# `prepare_intarna_out.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.intarna.prepare_predictions`.

Command: `coprarna helper prepare-intarna` followed by `srna up down region genbank-group`.

Regression fixtures: `tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
