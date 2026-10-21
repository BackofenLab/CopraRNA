# `parse_region_from_genome.pl` migration mapping

Status: **ported**.

Python counterpart: `coprarna.genome.parse_regions`.

Command: `coprarna helper genome-regions` followed by `genbank ntup ntdown region`.

Regression fixtures: `region-circular-5utr`, `region-linear-3utr`, `region-linear-5utr`, `region-linear-cds`, `region-circular-3utr`, `region-circular-cds`.

The original source is retained unchanged. See `AGENTS.md` and `README.2.md`
for validation scope and compatibility decisions.
