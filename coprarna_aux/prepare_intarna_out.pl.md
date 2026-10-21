# `prepare_intarna_out.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.intarna.prepare_predictions`.

Command: `coprarna helper prepare-intarna` followed by `srna up down region genbank-group`.

Regression fixtures: `tests/test_workflow.py::test_full_native_prediction_matches_legacy`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
