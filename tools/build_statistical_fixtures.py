"""Generate deterministic inputs for the original statistical scripts."""

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = (
    (ROOT / "tests/fixtures/cluster-predictions/input/NC_000001_opt.intarna.csv")
    .read_text()
    .splitlines()[0]
)
options = "sRNA FASTA:input_sRNA.fa\nnt upstream:200\nnt downstream:100\nregion:5utr\nrelative clustersize:0.5\ncore count:1\nwin size:150\nmax bp dist:100\nCopraRNA1:0\nverbose:0\nwebsrv:0\nnooi:0\ntop count:200\nroot:1\nenrich:0\nnoclean:1\ncons:0\nooifilt:0\nversion:CopraRNA 2.1.3\n"
rng = random.Random(714)
organisms = ["NC_000001", "NC_000002", "NC_000003", "NC_000004"]
annotations = []
clusters = []
for i in range(120):
    common = rng.random()
    cells = []
    for j, org in enumerate(organisms):
        p = round(0.0001 + 0.9998 * (0.3 * common + 0.7 * rng.random()), 7)
        tag = f"g{i:03d}_{j}"
        missing = i % 7 == 0 and j == 3
        cells.append(
            ""
            if missing
            else f"{tag}(gene{i}|-12|{p}|3|12|4|13|GeneID:{1000 + i * 4 + j})"
        )
        if not missing:
            row = (
                [
                    tag,
                    "ncrna_" + org.lower(),
                    "ACGT",
                    "ACGT",
                    "AC",
                    "GU",
                    "((&))",
                    "((&))",
                    "3",
                    "12",
                    "4",
                    "13",
                    "((&))",
                    "((&))",
                    "-12",
                ]
                + ["0"] * 19
                + ["-1", str(p), str(i + 1)]
            )
            clusters.append(";".join(row))
    annotations.append(",".join(["0", "0", *cells, f"product{i}", "", "0"]))
for name, script, args, outputs in [
    (
        "join-copra2",
        "join_pvals_coprarna2.R",
        ["NC_000001"],
        [
            "CopraRNA2_final_all_balanced.csv",
            "CopraRNA2_final_all_ooi.csv",
            "CopraRNA2_final_all_evo.csv",
        ],
    ),
    (
        "join-copra1",
        "join_pvals_coprarna1.R",
        ["opt_tags.clustered"],
        ["CopraRNA1_with_pvsample.csv"],
    ),
]:
    folder = ROOT / "tests/fixtures" / name
    (folder / "input").mkdir(parents=True, exist_ok=True)
    spec = dict(
        language="r",
        script="coprarna_aux/" + script,
        args=args,
        outputs=outputs,
        rng_seed=731,
    )
    (folder / "case.json").write_text(json.dumps(spec, indent=2) + "\n")
    (folder / "input/CopraRNA_option_file.txt").write_text(options)
    (folder / "input/zscore.weight").write_text(
        "".join(
            f"{org.lower()};{(0.1, 0.2, 0.3, 0.4)[j]}\n"
            for j, org in enumerate(organisms)
        )
    )
    (folder / "input/ncrna.fa").write_text(
        "".join(f">ncRNA_{org}\nACGTACGTACGT\n" for org in organisms)
    )
    (folder / "input/opt_tags.clustered").write_text(
        HEADER.replace("id1;", "d1;", 1)
        + ";clusternumber\n"
        + "\n".join(clusters)
        + "\n"
    )
    (folder / "input/CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv").write_text(
        ",".join(
            [
                "fdr",
                "p-value",
                *organisms,
                "Annotation",
                "Additional homologs",
                "Amount sampled",
            ]
        )
        + "\n"
        + "\n".join(annotations)
        + "\n"
    )
