#!/usr/bin/env python3
"""Regenerate the auditable per-source migration inventory and mapping documents."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from coprarna.helpers import HELPERS  # noqa: E402

PLANNED = {
    "CopraRNA2.pl": "coprarna.pipeline.run_prediction",
    "homology_intaRNA.pl": "coprarna.homology.run_homology",
    "prepare_intarna_out.pl": "coprarna.intarna.prepare_predictions",
    "combine_clusters.pl": "coprarna.combination.combine_clusters",
    "get_refseq_from_refid.pl": "coprarna.ncbi.download_genbank",
    "parallelize_target_alignments.pl": "coprarna.alignments.align_targets",
    "add_pval_to_csv_evdfit.R": "coprarna.energy.add_energy_pvalues",
    "join_pvals_coprarna1.R": "coprarna.combination.join_copra1",
    "join_pvals_coprarna2.R": "coprarna.combination.join_copra2",
    "copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R": "coprarna.conservation.compute_conservation",
    "script_R_plots_8.R": "coprarna.plots.plot_interaction_regions",
    "evo_heatmap.R": "coprarna.plots.plot_evolutionary_heatmap",
    "extract_functional_enriched.R": "coprarna.enrichment.extract_enriched_terms",
    "DAVIDWebService_CopraRNA.py": "coprarna.david.coprarna_report",
    "DAVIDWebService_IntaRNA_chartReport.py": "coprarna.david.intarna_report",
    "prepare_output_for_websrv_new.pl": "coprarna.webserver.prepare_webserver_output",
    "rasterize.js": "coprarna.webserver.rasterize_heatmap",
    "build_kegg2refseq.pl": "coprarna.maintenance.build_organism_mapping",
    "add_date_omit_incompatible.sh": "coprarna.maintenance.write_availability_header",
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
        else:
            entry.update(
                planned_target=PLANNED[basename], oracle_fixtures=sorted(fixtures)
            )
            body = f"Planned counterpart: `{PLANNED[basename]}`.\n\nThis module/function has **not been implemented**. The legacy source remains\nthe implementation for this stage. It is not called by the new Python package."
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
                else "No regression fixtures have been captured for this stage yet."
            )
            + "\n\nThe original source is retained unchanged. See `AGENTS.md` and `README.2.md`\nfor parity limitations and remaining migration work.\n"
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
