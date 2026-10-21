# `DAVIDWebService_CopraRNA.py` migration mapping

Status: **ported**.

Python counterpart: `coprarna.david.coprarna_report`.

Command: `coprarna helper david-copra` followed by `predictions count`.

Regression fixtures: `tests/test_adapters.py::test_david_requests_and_stdout_match_untouched_source`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
