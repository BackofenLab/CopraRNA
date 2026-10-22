#!/usr/bin/env python3
"""Build small, reviewable inputs; expected files come only from legacy scripts."""

from pathlib import Path
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures"
HEADER = "id1;id2;seq1;seq2;subseq1;subseq2;subseqDP;subseqDB;start1;end1;start2;end2;hybridDP;hybridDB;E;ED1;ED2;Pu1;Pu2;E_init;E_loops;E_dangleL;E_dangleR;E_endL;E_endR;seedStart1;seedEnd1;seedStart2;seedEnd2;seedE;seedED1;seedED2;seedPu1;seedPu2;E_norm"


def case(name, script, args=(), files=None, outputs=(), modules=(), unordered=()):
    directory = FIXTURES / name
    (directory / "input").mkdir(parents=True, exist_ok=True)
    for filename, data in (files or {}).items():
        (directory / "input" / filename).write_text(data)
    spec = dict(
        script="coprarna_aux/" + script,
        language="r" if script.endswith(".R") else "perl",
        args=list(args),
        outputs=list(outputs),
    )
    if name == "archive-readme":
        spec["legacy_parity"] = False
    if name == "sampled-1":
        spec["legacy_args"] = list(args)
        spec["args"] = list(args)[:1]
    if modules:
        spec["modules"] = list(modules)
    if unordered:
        spec["unordered"] = list(unordered)
    (directory / "case.json").write_text(json.dumps(spec, indent=2) + "\n")
    return directory / "input"


def row(tag, org="NC_000001", energy="-12", pval=None):
    values = (
        [
            tag,
            "ncRNA_" + org,
            "ACGT",
            "ACGU",
            "AC",
            "GU",
            "((&))",
            "((&))",
            "3",
            "8",
            "2",
            "7",
            "((&))",
            "((&))",
            energy,
        ]
        + ["0"] * 19
        + ["-1"]
    )
    if pval is not None:
        values.append(pval)
    return ";".join(values) + "\n"


def genome(path, circular=True, genes=True, cds=True, rrna=True):
    # Nonperiodic sequence makes strand and boundary offsets visible.
    import random

    rng = random.Random(9157)
    rec = SeqRecord(
        Seq("".join(rng.choices("ACGT", k=1100))),
        id="NC_000001.1",
        name="NC_000001",
        description="Synthetic regression genome",
    )
    rec.annotations = dict(
        molecule_type="DNA",
        topology="circular" if circular else "linear",
        organism="Synthetic organism",
    )
    locations = [
        SimpleLocation(0, 90, strand=1),
        SimpleLocation(350, 510, strand=1),
        SimpleLocation(530, 730, strand=-1),
        SimpleLocation(1010, 1100, strand=-1),
        CompoundLocation(
            [SimpleLocation(850, 890, strand=1), SimpleLocation(910, 960, strand=1)]
        ),
        CompoundLocation(
            [SimpleLocation(160, 190, strand=-1), SimpleLocation(110, 140, strand=-1)]
        ),
    ]
    for i, location in enumerate(locations):
        q = dict(
            locus_tag=[f"LT{i}"],
            gene=[f"g({i}),name"],
            db_xref=[f"GeneID:{100 + i}"],
            product=["hypothetical protein" if i == 1 else f"product {i}, alpha; beta"],
            note=[f"note {i}, fallback"],
            translation=["MELQK"],
        )
        if genes:
            rec.features.append(
                SeqFeature(
                    location,
                    type="gene",
                    qualifiers={
                        k: v
                        for k, v in q.items()
                        if k in ("locus_tag", "gene", "db_xref")
                    },
                )
            )
        if cds:
            rec.features.append(SeqFeature(location, type="CDS", qualifiers=q))
    if rrna:
        rec.features.append(
            SeqFeature(
                SimpleLocation(760, 800, strand=-1),
                type="rRNA",
                qualifiers={"product": ["16S ribosomal RNA"]},
            )
        )
    SeqIO.write(rec, path, "genbank")


def main():
    fasta = ">aaaa:LT0 first title\nAC GT\nAC\n>bbbb:LT1 second title\nA\tTT\n>cccc:LT2\nGCTA\n"
    case(
        "fasta-gene",
        "fasta2genefile.pl",
        ["input.fa"],
        {"input.fa": fasta},
        ["input.fa.gene", "input.fa.tit"],
    )
    case(
        "format-fasta",
        "format_fasta.pl",
        ["input.fa"],
        {"input.fa": fasta},
        modules=["Bio::SeqIO"],
    )
    case(
        "blast-homology",
        "blast2homfile.pl",
        ["blast.tsv"],
        {
            "blast.tsv": "# comment\nb a 75 30 2 0 3 32 5 34 1e-9 20\na b 80 30 2 0 2 31 8 37 1e-8 30\na c 10 30 2 0 1 30 1 30 1e-7 40\nc c 100 30 0 0 1 30 1 30 0 50\na d 90 30 1 0 1 30 1 30 0.1 15\n"
        },
    )
    case(
        "transform-distance",
        "transform_distmat.pl",
        ["distmat.out"],
        {
            "distmat.out": "header\n" * 8
            + "\t0.00\t12.30\t4.56\tNC_000001 1\n\t0.00\t7.89\tNC_000002 2\n\t0.00\tNZ_ABC3 3\n"
        },
    )
    for zero in (False, True):
        tree = (
            "((nc_000001:0.10000,nc_000002:0.10000):0.20000,nz_abc3:0.30000);\n"
            if not zero
            else "(nc_000001:0.00000,nc_000002:0.00000);\n"
        )
        branches = (
            "2 nc_000001 0.10000 0.0\n2 nc_000002 0.10000 0.0\n1 2 -0.20000 0.1\n1 nz_abc3 0.30000 0.0\n"
            if not zero
            else "1 nc_000001 0.00000 0.0\n1 nc_000002 0.00000 0.0\n"
        )
        case(
            "weights-zero" if zero else "weights-branches",
            "compute_weights.pl",
            ["tree", "neighbor"],
            {
                "tree": tree,
                "neighbor": f"{2 if zero else 3} Populations\nFrom To Length Height\n"
                + branches,
            },
        )
    for cop2 in ("1",):
        case(
            "sampled-" + cop2,
            "get_amount_sampled_values_and_add_to_table.pl",
            ["result.csv", cop2],
            {
                "result.csv": "fdr,p-value,org1,org2,Annotation,Additional homologs\n0.1,0.01,tag,,note,\n0.2,0.02,,tag,,\n"
            },
        )
    enrichment = 'score=2.5\nafdr=0.05\ncategoryName="GO"\ntermName="a term"\nlistHits=2\npercent=10\nease=0.01\ngeneIds="100, 101"\nlistTotals=20\npopHits=4\npopTotals=500\nfoldEnrichment=3\nbonferroni=0.2\nbenjamini=0.1\n'
    case(
        "enrichment-report",
        "make_enrichment_table_from_py_output.pl",
        ["report.txt"],
        {"report.txt": enrichment},
    )
    case(
        "additional-homologs",
        "parse_homologs_from_domclust_table.pl",
        ["result.csv", "cluster.tab"],
        {
            "result.csv": "p-value,NC_000001,Annotation\n0.1,lt0(g|GeneID:100),note\n0.2,,none\n",
            "cluster.tab": "# header\n1\taaaa:LT0(1) aaaa:LT1(2)\tbbbb:LT2(1)\n",
        },
        modules=["List::MoreUtils"],
    )
    option = "sRNA FASTA:input_sRNA.fa\nnt upstream:200\nnt downstream:100\nregion:5utr\nrelative clustersize:0.5\n"
    data = (
        HEADER
        + ";p-value\n"
        + row("LT0", pval="0.2")
        + row("LT1", pval="0.1")
        + row("LT2", "NC_000002", pval="0.8")
    )
    cluster = "#cluster\taaaa\tbbbb\n1\taaaa:LT0(1) aaaa:LT1(2)\tbbbb:LT2(1)\n2\taaaa:missing\tbbbb:LT0\n"
    case(
        "cluster-predictions",
        "cluster_intarna_csv.pl",
        files={"NC_000001_opt.intarna.csv": data, "cluster.tab": cluster},
        modules=["List::MoreUtils"],
    )
    for mode in ("circular", "linear"):
        for region in ("5utr", "3utr", "cds"):
            folder = case(
                f"region-{mode}-{region}",
                "parse_region_from_genome.pl",
                ["NC_000001.gb", "200", "100", region],
                modules=["Bio::SeqIO"],
            )
            genome(folder / "NC_000001.gb", circular=mode == "circular")
    for genes, cds in ((True, True), (True, False), (False, False)):
        folder = case(
            f"check-features-{int(genes)}-{int(cds)}",
            "check_for_gene_CDS_features.pl",
            ["NC_000001.gb"],
            modules=["Bio::SeqIO"],
        )
        genome(folder / "NC_000001.gb", genes=genes, cds=cds)
    folder = case(
        "extract-proteins",
        "get_CDS_from_gbk.pl",
        ["NC_000913.gb"],
        modules=["Bio::SeqIO"],
    )
    genome(folder / "NC_000913.gb")
    folder = case(
        "extract-16s",
        "parse_16s_from_gbk.pl",
        ["NC_000001.gb,NC_000002.gb"],
        modules=["Bio::SeqIO"],
    )
    genome(folder / "NC_000001.gb", rrna=False)
    genome(folder / "NC_000002.gb")
    folder = case(
        "annotate-intarna",
        "add_GI_genename_annotation_intarna.pl",
        files={"intarna_websrv_table.csv": data},
        outputs=["intarna_websrv_table_ncbi.csv"],
        modules=["Bio::SeqIO"],
    )
    genome(folder / "NC_000001.gb")
    opt = "NC_000001_opt.intarna.csv"
    sub = "NC_000001_subopt.intarna.csv"
    case(
        "disentangle",
        "disentangle_sorted_intarna_CSV.pl",
        files={
            "NC_000001.fa.intarna.sorted.csv": HEADER
            + "\n"
            + row("LT0")
            + row("LT1")
            + row("LT0", energy="-10")
            + row("LT0", energy="-8")
        },
        outputs=[opt, sub],
        unordered=[opt, sub],
    )
    case(
        "adjust-pvalues",
        "calc_padj.R",
        ["pvalues"],
        {"pvalues": "p\n0\n0.01\n0.1\n0.1\n1\nNA\n"},
        ["padj.csv"],
    )
    clustered = (
        HEADER.replace("id1;", "d1;", 1)
        + ";p-value;clusternumber\n"
        + row("lt0", "NC_000002", pval="0.2").rstrip()
        + ";1\n"
        + row("lt1", "NC_000001", pval="0.1").rstrip()
        + ";1\n"
        + row("lt2", "NC_000002", pval="0.8").rstrip()
        + ";2\n"
    )
    case(
        "prepare-copra2",
        "prep_cop2.R",
        ["opt_tags.clustered"],
        {"opt_tags.clustered": clustered},
        ["CopraRNA2_prep.csv"],
    )
    folder = case(
        "annotate-raw",
        "annotate_raw_output.pl",
        ["result.raw", "opt_tags.clustered", "NC_000001.gb", "NC_000002.gb"],
        files={
            "result.raw": "0.001;other0;lt0;\n0.2;lt1;\n0.3;other1;\n",
            "opt_tags.clustered": clustered
            + row("other0", "NC_000002", pval="0.05").rstrip()
            + ";1\n"
            + row("other1", "NC_000002", pval="0.15").rstrip()
            + ";2\n",
        },
        modules=["Bio::SeqIO"],
    )
    genome(folder / "NC_000001.gb")
    genome(folder / "NC_000002.gb")
    text = (folder / "NC_000002.gb").read_text()
    for i in range(6):
        text = text.replace(f"LT{i}", f"other{i}")
    (folder / "NC_000002.gb").write_text(text)
    case(
        "heatmap-json",
        "make_heatmap_json.pl",
        ["enrichment.tsv"],
        {
            "enrichment.tsv": "2.5\t1.2\n2\t2\t1\t2\n"
            + "\t".join(
                [
                    "x0",
                    "p",
                    "x2",
                    "x3",
                    "source",
                    "x5",
                    "x6",
                    "x7",
                    "x8",
                    "x9",
                    "x10",
                    "geneA",
                    "geneB",
                ]
            )
            + "\n"
            + "\t".join(
                [
                    "a",
                    "0.001",
                    "b",
                    "c",
                    "GO_term",
                    "d",
                    "e",
                    "f",
                    "g",
                    "h",
                    "i",
                    "1",
                    "0",
                ]
            )
            + "\n"
            + "\t".join(
                [
                    "a",
                    "0.1",
                    "b",
                    "c",
                    "Other_term",
                    "d",
                    "e",
                    "f",
                    "g",
                    "h",
                    "i",
                    "0",
                    "1",
                ]
            )
            + "\n"
        },
        outputs=["copraRNA.json"],
    )
    case(
        "specific-targets",
        "find_single_specific_targets_in_termCluster.pl",
        files={
            "termClusterReport.txt": "Annotation Cluster 1\tEnrichment Score:  2.5\nCategory\tTerm\tCount\t%\tPvalue\tGenes\nGO\ttermA\t1\t10\t0.01\t100\n",
            "IntaRNA_chartReport_grepped.txt": "geneIds=100,101\ntermName=termA\n",
            "intarna_websrv_table_ncbi.csv": HEADER
            + ";p-value;genename;entrez_gene_id;annotation\n"
            + row("LT1", pval="0.1").rstrip()
            + ";gene1;101;annotation\n",
        },
    )
    (case("archive-readme", "print_archive_README.pl") / ".gitkeep").touch()
    case(
        "adjust-pvalues-precision",
        "calc_padj.R",
        ["pvalues"],
        {
            "pvalues": "p\n0.0001\n0.00015\n1e-15\n1e-5\n0.00000012345\n0.1234567891234567\n1\n0\n"
        },
        ["padj.csv"],
    )
    folder = case(
        "fit-energy-pvalues",
        "add_pval_to_csv_evdfit.R",
        files={
            "CopraRNA_option_file.txt": option,
            opt: HEADER
            + "\n"
            + "".join(
                row("lt" + str(i), energy=str(-5 - i * 0.31 - (i % 7)))
                for i in range(80)
            ),
            sub: HEADER
            + "\n"
            + "".join(
                row("lt" + str(i), energy=str(-3 - i * 0.25 - (i % 5)))
                for i in range(80)
            ),
        },
        outputs=[opt, sub],
    )


if __name__ == "__main__":
    main()
