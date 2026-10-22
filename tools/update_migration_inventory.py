#!/usr/bin/env python3
"""Regenerate the auditable per-source migration inventory and mapping documents."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from coprarna.helpers import HELPERS  # noqa: E402

RETIRED = {"join_pvals_coprarna1.R", "remove_clusters_under_percantage.R"}

COVERAGE = {
    "CopraRNA2.pl": "tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors",
    "homology_intaRNA.pl": "tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors",
    "prepare_intarna_out.pl": "tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors",
    "combine_clusters.pl": "tests/test_workflow.py::test_corrected_workflow_builds_clusters_preserves_inputs_and_exports_vectors",
    "parallelize_target_alignments.pl": "tests/test_workflow.py::test_parallel_alignment_inputs_and_outputs_match",
    "copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R": "tests/test_workflow.py::test_conservation_and_all_203_jalview_exports_match",
    "script_R_plots_8.R": "tests/test_workflow.py::test_plot_geometry_and_heatmap_match_original_R",
    "evo_heatmap.R": "tests/test_workflow.py::test_plot_geometry_and_heatmap_match_original_R",
    "DAVIDWebService_CopraRNA.py": "tests/test_adapters.py::test_david_requests_and_stdout_match_untouched_source",
    "DAVIDWebService_IntaRNA_chartReport.py": "tests/test_adapters.py::test_david_requests_and_stdout_match_untouched_source",
    "rasterize.js": "tests/test_adapters.py::test_rasterize_real_chromium",
    "build_kegg2refseq.pl": "tests/test_adapters.py::test_maintenance_mapping_matches_seeded_legacy",
    "add_date_omit_incompatible.sh": "tests/test_adapters.py::test_maintenance_mapping_matches_seeded_legacy",
}


def main():
    original = json.loads((ROOT / "tests/fixtures/legacy-sha256.json").read_text())
    cases = {
        p.parent.name: json.loads(p.read_text())
        for p in (ROOT / "tests/fixtures").glob("*/case.json")
    }
    by_legacy = {helper.legacy: (name, helper) for name, helper in HELPERS.items()}
    inventory = []
    for source in sorted(original):
        if Path(source).suffix not in (".pl", ".R", ".py", ".sh", ".js"):
            continue
        basename = Path(source).name
        entry = {"legacy": source, "status": "pending"}
        fixtures = [name for name, spec in cases.items() if spec["script"] == source]
        if basename in COVERAGE:
            fixtures.append(COVERAGE[basename])
        if basename in by_legacy:
            name, helper = by_legacy[basename]
            entry.update(
                status="ported",
                target=helper.target,
                command="coprarna helper " + name,
                fixtures=sorted(fixtures),
            )
            body = f"Python counterpart: `{helper.target}`.\n\nCommand: `coprarna helper {name}`"
            if helper.arguments:
                body += " followed by `" + " ".join(helper.arguments) + "`"
            body += "."
        elif basename in RETIRED:
            entry.update(
                status="retired",
                reason="Not part of the CopraRNA2-only runtime requested in review",
            )
            body = "Retired from the Python implementation under the CopraRNA2-only scope. No runtime counterpart is provided."
        else:
            raise ValueError(f"Unmapped original source: {source}")
        mapping = ROOT / (source + ".md")
        mapping.write_text(
            f"# `{basename}` migration mapping\n\nStatus: **{entry['status']}**.\n\n"
            + body
            + "\n\n"
            + (
                "Regression fixtures: "
                + ", ".join(f"`{name}`" for name in fixtures)
                + "."
                if fixtures
                else "Retired from the current runtime; no regression fixture is retained."
            )
            + "\n\nThe original source is retained unchanged. See `AGENTS.md` and `README.2.md`\nfor validation scope and compatibility decisions.\n"
        )
        inventory.append(entry)
    (ROOT / "coprarna/migration.json").write_text(
        json.dumps(inventory, indent=2) + "\n"
    )
    print(
        f"{sum(item['status'] == 'ported' for item in inventory)}/{len(inventory)} scripts ported"
    )


if __name__ == "__main__":
    main()
