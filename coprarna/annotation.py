"""GenBank annotation helpers with the original field ordering and sanitization."""

from pathlib import Path
import re

from .io import genbank_records, genome_paths, accession


def _annotation(qualifiers):
    """Prefer a non-hypothetical product, then note/function annotation."""
    product = qualifiers.get("product", [""])[0]
    if "product" in qualifiers and "hypothetical protein" not in product:
        return product
    return qualifiers.get("note", qualifiers.get("function", [""]))[0]


def annotate_intarna(directory: Path) -> str:
    """Add gene names, GeneIDs and sanitized annotations from compressed genomes."""
    lookup = {}
    for path in genome_paths(directory):
        for record in genbank_records(path):
            for feature in record.features:
                q = feature.qualifiers
                if feature.type != "CDS" or "locus_tag" not in q:
                    continue
                gid = ""
                for value in q.get("db_xref", []):
                    if "GeneID:" in value:
                        gid = value.replace("GeneID:", "")
                annotation = _annotation(q).replace(",", "").replace(";", "")
                lookup[q["locus_tag"][0]] = ";".join(
                    [q.get("gene", [""])[0], gid, annotation]
                )
    output = []
    for line in (directory / "intarna_websrv_table.csv").read_text().splitlines():
        if re.search(
            "id1;id2;seq1;seq2;subseq1;subseq2;subseqDP;subseqDB;start1;end1;start2;end2;hybridDP;",
            line,
        ):
            output.append(line + ";genename;entrez_gene_id;annotation\n")
        else:
            output.append(line + ";" + lookup.get(line.split(";")[0], ";;") + "\n")
    (directory / "intarna_websrv_table_ncbi.csv").write_text("".join(output))
    return ""


def annotate_raw(
    result: Path, clustered: Path, groups: list[str], directory: Path
) -> str:
    """Expand cluster tags into organism columns with interaction and gene metadata.

    Groups specify each organism's replicons, primary first. Keep the historical
    field order and text sanitization because downstream tables parse these cells."""
    names, gids, columns, products = {}, {}, {}, {}
    header = ["p-value"]
    for col, group in enumerate(groups, 2):
        header.append(accession(group.split(",")[0]))
        for filename in group.split(","):
            record = next(genbank_records(directory / filename))
            for feature in record.features:
                q = feature.qualifiers
                ids = list(q.get("locus_tag", []))
                if not ids:
                    ids.append("N/A")
                ids += q.get("gene", [])
                if len(ids) < 2:
                    ids.append("N/A")
                ids += q.get("db_xref", [])
                if len(ids) < 3:
                    ids.append("N/A")
                tag, name = ids[0].lower(), ids[1]
                gid = next(
                    (value for value in reversed(ids[2:]) if "GeneID:" in value), ""
                )
                # The lowercase N/A comparison is retained from Perl.
                if tag != "N/A":
                    columns[tag] = col
                if tag != "N/A" and name != "N/A":
                    names[tag] = (
                        name.replace(",", ";").replace("(", "").replace(")", "")
                    )
                if tag != "N/A" and gid != "N/A" and "GeneID:" in gid:
                    gids[tag] = gid
    for filename in groups[0].split(","):
        for record in genbank_records(directory / filename):
            for feature in record.features:
                q = feature.qualifiers
                if feature.type == "CDS" and "locus_tag" in q:
                    value = _annotation(q)
                    if value:
                        products[q["locus_tag"][0]] = value
    properties = {}
    for line in clustered.read_text().splitlines():
        row = line.split(";")
        properties[row[0].lower()] = "|" + "|".join(
            row[i] for i in (14, -2, 8, 9, 10, 11)
        )
    output = [",".join(header) + ",Annotation\n"]
    for line in result.read_text().splitlines():
        cells = line.split(";")[:-1]
        by_column = {columns.get(tag, 0): tag for tag in cells[1:]}
        tags = [by_column[index] for index in sorted(by_column)]
        row = [cells[0]] + [""] * len(groups)
        for col, tag in by_column.items():
            if 2 <= col <= len(groups) + 1:
                row[col - 1] = (
                    f"{tag}({names.get(tag, 'N/A')}{properties.get(tag, '')}|{gids.get(tag, 'N/A')})"
                )
        product = next(
            (
                value
                for key, value in products.items()
                if tags and re.search(tags[0], key, re.I)
            ),
            "",
        )
        output.append(",".join(row) + "," + product.replace(",", " ") + "\n")
    return "".join(output)
