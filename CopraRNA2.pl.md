# `CopraRNA2.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.pipeline.run_prediction`.

Command: `coprarna helper predict` followed by `srna`.

Regression fixtures: `tests/test_workflow.py::test_full_native_prediction_matches_legacy`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
