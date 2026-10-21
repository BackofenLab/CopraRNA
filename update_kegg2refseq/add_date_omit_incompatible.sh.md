# `add_date_omit_incompatible.sh` migration mapping

Status: **ported**.

Python counterpart: `coprarna.maintenance.write_availability_header`.

Command: `coprarna helper availability-header` followed by `header-file`.

Regression fixtures: `tests/test_adapters.py::test_maintenance_mapping_matches_seeded_legacy`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
