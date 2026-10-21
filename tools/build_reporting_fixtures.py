"""Build deterministic inputs; expected outputs are captured from original scripts."""

from pathlib import Path
import json
from tests.support.legacy import FIXTURES, capture


def make(name, script, args, outputs, language="perl"):
    case = FIXTURES / name
    case.mkdir(exist_ok=True)
    (case / "input").mkdir(exist_ok=True)
    (case / "expected").mkdir(exist_ok=True)
    (case / "case.json").write_text(
        json.dumps(
            dict(language=language, script=script, args=args, outputs=outputs), indent=2
        )
        + "\n"
    )
    return case


def build():
    case = make(
        "functional-enrichment",
        "coprarna_aux/extract_functional_enriched.R",
        ["prediction.csv", "report.txt", "enrichment.txt"],
        ["enrichment.txt"],
        "r",
    )
    header = "fdr,p-value,NC_000913,Annotation,Additional homologs,Amount sampled\n"
    rows = [
        f".01,.002,tag{i}(gene{i}|-12|0.01|2|20|3|15|GeneID:{100 + i}),test,,0\n"
        for i in range(300)
    ]
    (case / "input/prediction.csv").write_text(header + "".join(rows))
    (case / "input/report.txt").write_bytes(
        (FIXTURES / "enrichment-report/expected/stdout").read_bytes()
    )
    case = make(
        "webserver-output",
        "coprarna_aux/prepare_output_for_websrv_new.pl",
        ["prediction.csv", "interactions.txt"],
        ["coprarna_internal_table.csv"],
    )
    (case / "input/prediction.csv").write_text(header + "".join(rows[:2]))
    (case / "input/CopraRNA_option_file.txt").write_text("CopraRNA1:0\n")
    ui = []
    for i in range(2):
        ui.extend(
            [
                f">tag{i}",
                *[f"header{k}" for k in range(4)],
                "5'-AAA       UUU-3'",
                "ACGUACGUACGUACGUACG",
                "UGCAUGCAUGCAUGCAUGC",
                "3'-CCC       GGG-5'",
                "",
                "seed target: 3 -- 9",
                "",
                "",
                "seed query: 5 -- 11",
                "",
                "unfolding target: 1.2",
                "unfolding query: 2.3",
                "hybridization: -15.5",
            ]
        )
    (case / "input/interactions.txt").write_text("\n".join(ui) + "\n")
    case = make(
        "cached-genbank",
        "coprarna_aux/get_refseq_from_refid.pl",
        ["-acc", "NC_000913", "-g", "cached.gb"],
        ["cached.gb"],
    )
    spec = json.loads((case / "case.json").read_text())
    spec["modules"] = ["Bio::DB::EUtilities"]
    (case / "case.json").write_text(json.dumps(spec, indent=2) + "\n")
    (case / "input/cached.gb").write_text("cached genome bytes\n")
    import tempfile

    for name in ("functional-enrichment", "webserver-output", "cached-genbank"):
        case = FIXTURES / name
        with tempfile.TemporaryDirectory() as tmp:
            for name, data in capture(case, Path(tmp)).items():
                dest = case / "expected" / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)


if __name__ == "__main__":
    build()
