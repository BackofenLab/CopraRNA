# `rasterize.js` migration mapping

Status: **ported**.

Python counterpart: `coprarna.webserver.rasterize_heatmap`.

Command: `coprarna helper rasterize` followed by `address output`.

Regression fixtures: `tests/test_adapters.py::test_rasterize_real_chromium`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
