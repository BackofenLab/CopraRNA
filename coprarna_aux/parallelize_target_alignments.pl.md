# `parallelize_target_alignments.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.alignments.align_targets`.

Command: `coprarna helper align-targets` followed by `prediction-table`.

Regression fixtures: `tests/test_workflow.py::test_parallel_alignment_inputs_and_outputs_match`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
