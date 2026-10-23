"""Webserver table formatting and optional browser report exports."""

from pathlib import Path
import re
import json
from .runtime import command
from . import parameters

HEADER = "Rank,CopraRNA p-value,CopraRNA fdr,Locus Tag,Gene Name,Energy kcal/mol,IntaRNA p-value,Position mRNA,Position ncRNA,Annotation,Additional homologs,Entrez GeneID,Interaction,Position Seed - mRNA,Position Seed - ncRNA,Hybridization Energy kcal/mol,Unfolding Energy - mRNA kcal/mol,Unfolding Energy - ncRNA kcal/mol,Amount sampled"


def prune_interaction(line1, line2, line3, line4, target, query):
    """Format a compact four-line duplex into the webserver's escaped text field."""
    line1, line2, line3, line4 = (
        value.strip() for value in (line1, line2, line3, line4)
    )
    left, *_, right = re.split(r"\s", line1)
    bottomleft, *_, bottomright = re.split(r"\s", line4)
    tstart, tend = map(int, re.findall(r"\d+", target))
    qstart, qend = map(int, re.findall(r"\d+", query))
    topflags = int(len(left) == 3) + int(len(right) == 3)
    tl = " " * (14 if len(left) == 3 else 13)
    pstart = tstart - 1 + int(len(left) == 3)
    pend = tend + 1 - int(len(right) == 3)
    output = (
        tl
        + str(pstart)
        + " " * max(0, 1 + len(line2) - len(str(pstart)) - topflags)
        + str(pend)
        + "\\n"
    )
    output += tl + "|" + " " * max(0, len(line2) - topflags) + "|\\n"

    def shorten(value, isleft):
        if len(value) > 11:
            return (
                value[:6] + "..." + value[-5:]
                if isleft
                else value[:5] + "..." + value[-6:]
            )
        return value.rjust(14) if isleft else value

    topleft = shorten(left, True)
    middle = re.sub(right, "", re.sub(left, "", line1, count=1), count=1)
    output += topleft + middle + shorten(right, False) + "\\n"
    output += " " * len(topleft) + line2 + "\\n"
    bonds = "".join(
        "|"
        if a + b in ("AU", "UA", "GC", "CG")
        else ":"
        if a + b in ("GU", "UG")
        else " "
        for a, b in zip(line2, line3)
    )
    output += " " * len(topleft) + bonds + "\\n" + " " * len(topleft) + line3 + "\\n"
    middle = re.sub(bottomright, "", re.sub(bottomleft, "", line4, count=1), count=1)
    output += shorten(bottomleft, True) + middle + shorten(bottomright, False) + "\\n"
    bottomflags = int(len(bottomleft) == 3) + int(len(bottomright) == 3)
    bl = " " * (14 if len(bottomleft) == 3 else 13)
    output += bl + "|" + " " * max(0, len(line3) - bottomflags) + "|\\n"
    qstart = qstart - 1 + int(len(bottomright) == 3)
    qend = qend + 1 - int(len(bottomleft) == 3)
    output += (
        bl
        + str(qend)
        + " " * max(0, 1 + len(line3) - len(str(qend)) - bottomflags)
        + str(qstart)
        + "\\n"
    )
    return re.sub(r"\s", "_", output)


def prepare_webserver_output(result: Path, interactions: Path, directory: Path) -> str:
    """Write the frontend table from comparative rows and native-derived duplexes.

    Keep rows with missing comparative probabilities as NA; the more selective
    sTar personality can leave valid organism interactions without a combined score."""
    lines = result.read_text().splitlines()
    ui = interactions.read_text().splitlines()
    if len(lines) <= 1:
        raise ValueError(f"No data in {result}")
    output = [HEADER + "\n"]
    rank = 0
    for line in lines[1:]:
        row = line.split(",")
        if not row[2]:
            continue
        rank += 1
        cell = re.split(r"[\(,|]", row[2])
        tag = cell[0]
        target = cell[4] + " -- " + cell[5]
        query = cell[6] + " -- " + cell[7]
        match = re.search(r"GeneID:(\d+)\)", cell[8])
        gid = match[1] if match else ""
        start = next(
            (i for i, line in enumerate(ui) if line.casefold() == ">" + tag.casefold()),
            None,
        )
        block = [""] + (ui[start + 5 : start + 19] if start is not None else [""] * 14)
        block += [""] * (15 - len(block))
        for i in (7, 10, 12, 13, 14):
            block[i] = re.sub(".*: ", "", block[i])
        interaction = prune_interaction(*block[1:5], target, query)

        def number(value, precision):
            return "NA" if value in ("", "NA") else format(float(value), precision)

        fields = [
            str(rank),
            number(row[1], ".4g"),
            number(row[0], ".4g"),
            tag,
            cell[1],
            number(cell[2], ".2f"),
            number(cell[3], ".6f"),
            target,
            query,
            row[-3],
            row[-2],
            gid,
            interaction,
            block[7],
            block[10],
            block[14],
            block[12],
            block[13],
            "0",
        ]
        output.append(",".join(fields) + "\n")
    (directory / "coprarna_internal_table.csv").write_text("".join(output))
    return ""


def format_interaction_ui(row):
    """IntaRNA 2.2 OutputHandlerIntaRNA1 layout, using current native CSV data.

    Reference: BackofenLab/IntaRNA v2.2.0, OutputHandlerIntaRNA1.cpp.
    The removed Perl UI wrapper supplied only these scientific command options.
    """
    tag, target, query_id, query, pairing, energy, ed1, ed2, *seeds = row
    pairs = [
        (int(a) - 1, int(b) - 1) for a, b in re.findall(r"\((\d+),(\d+)\)", pairing)
    ]
    if not pairs:
        return ""
    pairs.sort()
    left, right = pairs[0], pairs[-1]
    flank = 3 + max(left[0], len(query) - left[1] - 1)
    top = ("5'-" + target[: left[0]]).rjust(flank)
    bottom = ("3'-" + query[left[1] + 1 :][::-1]).rjust(flank)
    upper = " " * flank
    lower = " " * flank
    for k, (t, q) in enumerate(pairs):
        if k:
            previous_t, previous_q = pairs[k - 1]
            one = target[previous_t + 1 : t]
            two = query[q + 1 : previous_q][::-1]
            width = max(len(one), len(two))
            top += one.ljust(width)
            bottom += two.ljust(width)
            upper += " " * width
            lower += " " * width
        top += " "
        bottom += " "
        upper += target[t]
        lower += query[q]
    top += target[right[0] + 1 :] + "-3'"
    bottom += query[: right[1]][::-1] + "-5'"
    seed = [value.split(":")[0] for value in seeds]

    def fmt(value):
        return format(float(value), ".6g")

    return (
        f">{tag}\n{target}\n>{query_id}\n{query}\n\n"
        + top
        + "\n"
        + upper
        + "\n"
        + lower
        + "\n"
        + bottom
        + "\n\n"
        + f"positions(target)     : {left[0] + 1} -- {right[0] + 1}\npositions seed(target): {seed[0]} -- {seed[1]}\npositions with dangle(target): {left[0] + 1} -- {right[0] + 1}\npositions(ncRNA)      : {right[1] + 1} -- {left[1] + 1}\npositions seed(ncRNA) : {seed[2]} -- {seed[3]}\npositions with dangle(ncRNA): {right[1] + 1} -- {left[1] + 1}\nED target need: {fmt(ed1)} kcal/mol\nED ncRNA  need: {fmt(ed2)} kcal/mol\nhybrid energy : {fmt(float(energy) - float(ed1) - float(ed2))} kcal/mol\n\nenergy: {fmt(energy)} kcal/mol\n"
    )


def export_prediction(directory, ooi, up, down):
    """Predict selected OOI targets as native gzip CSV and build frontend exports.

    Forward the same parameter file and interaction count as the primary run.
    Apply the same artificial-helix filter before choosing one retained native
    interaction per target for the frontend. This prevents a removed best hit
    from reappearing when the comparative result used a retained suboptimal hit.
    Require the main run's context sidecar, reuse its extended genomic sequences
    and restrict/rebase interactions to the requested region. A missing or stale
    sidecar fails before native prediction rather than losing accessibility context.
    """
    from Bio import SeqIO
    from .genome import TargetRegion
    from .intarna import _predict_targets

    parameters.check_version(directory)
    path = directory / f"{ooi}_upfromstartpos_{up}_down_{down}.fa"
    records = list(SeqIO.parse(path, "fasta"))
    lines = (directory / "CopraRNA_result.csv").read_text().splitlines()[1:]
    tags = [line.split(",")[2].split("(")[0] for line in lines if line.split(",")[2]]
    context_path = Path(str(path) + ".context.json")
    if not context_path.is_file():
        raise ValueError("Missing genomic accessibility context; rerun prepare-intarna")
    context = {
        r["tag"]: TargetRegion(**r) for r in json.loads(context_path.read_text())
    }
    selected = []
    for tag in dict.fromkeys(t.casefold() for t in tags):
        for record in records:
            if tag == record.id.casefold():
                target = context[record.id]
                if target.requested != str(record.seq):
                    raise ValueError(
                        f"Stale genomic accessibility context for {record.id}"
                    )
                selected.append(target)
    (directory / "CopraRNA_top_targets.fa").write_text(
        "".join(
            f">{r.id}\n{r.seq}\n"
            for tag in tags
            for r in records
            if tag.casefold() == r.id.casefold()
        )
    )
    columns = "id1,seq1,id2,seq2,bpList,E,ED1,ED2,seedStart1,seedEnd1,seedStart2,seedEnd2,hybridDP,start1,end1,start2,end2"
    _predict_targets(
        selected,
        f"ncRNA_{ooi}.fa",
        "IntaRNA_web.csv.gz",
        directory,
        columns,
        runner=command,
    )
    from .rcompat import read_table

    from .intarna import filter_perfect_helices

    output = directory / "IntaRNA_web.csv.gz"
    audit = filter_perfect_helices(output, parameters.values(directory))
    (directory / "CopraRNA_web.filter.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    _, rows = read_table(output, sep=";")
    retained = {}
    for row in rows:
        tag = row[0].casefold()
        if tag not in retained or float(row[5]) < float(retained[tag][5]):
            retained[tag] = row
    (directory / "IntaRNA_web.txt").write_text(
        "".join(format_interaction_ui(row[:12]) for row in retained.values())
    )
    prepare_webserver_output(
        directory / "CopraRNA_result.csv",
        directory / "IntaRNA_web.txt",
        directory,
    )
    (directory / "coprarna_internal_table.csv").rename(
        directory / "coprarna_websrv_table.csv"
    )
    import shutil

    shutil.copyfile(path, directory / "target_sequences_orgofint.fa")


def rasterize_heatmap(address: str, output: Path, paper=None, zoom=1.0) -> str:
    """Optionally export a browser report as PDF, SVG or a raster image.

    Prediction uses Python vector rendering directly; this browser adapter is for
    explicit helper use and requires the optional Playwright/Chromium installation."""
    from playwright.sync_api import sync_playwright

    with sync_playwright() as browser_runtime:
        browser = browser_runtime.chromium.launch(headless=True)
        page = browser.new_page(
            viewport={"width": 600, "height": 600}, device_scale_factor=float(zoom)
        )
        url = address if "://" in address else Path(address).resolve().as_uri()
        page.goto(url, wait_until="networkidle")
        page.wait_for_timeout(200)
        if output.suffix == ".svg":
            from .plots import write_svg_pages

            drawings = []
            for frame in page.frames:
                drawings.extend(
                    frame.locator("svg").evaluate_all(
                        "(nodes) => nodes.map(node => new XMLSerializer().serializeToString(node))"
                    )
                )
            if not drawings:
                raise ValueError("The report contains no SVG graphics")
            write_svg_pages(drawings, output)
        elif output.suffix == ".pdf":
            settings = {"format": paper or "A4", "print_background": True}
            if paper and "*" in paper:
                width, height = paper.split("*")
                settings = dict(width=width, height=height, print_background=True)
            page.pdf(path=str(output), **settings)
        else:
            page.screenshot(path=str(output))
        browser.close()
    return ""
