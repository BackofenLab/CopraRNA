# `CopraRNA2.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.pipeline.run_prediction`.

Command: `coprarna helper predict` followed by `srna`.

Regression fixtures: `tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
