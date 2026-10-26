"""GenBank helpers preserving the legacy coordinate and boundary conventions."""

from pathlib import Path
from dataclasses import dataclass
import re

from .io import genbank_records, accession


def _reverse_complement(sequence: str) -> str:
    from Bio.Seq import Seq

    return str(Seq(sequence).reverse_complement())


def _subseq(sequence: str, start: int, end: int) -> str:
    # Bio::Seq allows an empty interval when start == end+1.
    if start < 1 or end > len(sequence) or start > end + 1:
        raise ValueError(f"Invalid subsequence coordinates: {start}, {end}")
    return sequence[start - 1 : end]


def get_region(
    sequence: str,
    circular: bool,
    strand: int,
    pos: int,
    upstream: int,
    downstream: int,
    tag: str,
    *,
    cds: bool = False,
) -> str:
    """Direct port of get_region/get_region_cds; positions are 1-based."""
    length = len(sequence)

    def sub(start, end):
        return _subseq(sequence, start, end)

    truncated = False
    if strand == -1:
        if pos - downstream >= 1 and pos + upstream <= length:
            value = _reverse_complement(sub(pos - downstream + 1, pos + upstream))
        elif not circular:
            truncated = True
            if pos - downstream < 1:
                value = _reverse_complement(sub(1, pos + upstream))
            else:
                value = _reverse_complement(sub(pos - downstream + 1, length))
        elif pos - downstream < 1:
            overhang = downstream - pos
            value = _reverse_complement(sub(1, pos + upstream)) + _reverse_complement(
                sub(length - overhang + 1, length)
            )
        else:
            overhang = pos + upstream - length
            value = _reverse_complement(sub(1, overhang)) + _reverse_complement(
                sub(pos - downstream + 1, length)
            )
    else:
        if pos - upstream >= 1 and pos + downstream <= length:
            value = sub(pos - upstream, pos + downstream - 1)
        elif not circular:
            truncated = True
            if pos - upstream < 1:
                value = sub(1, pos + downstream - 1)
            else:
                value = sub(pos - upstream, length)
        elif pos - upstream < 1:
            overhang = upstream - pos
            value = sub(length - overhang, length) + sub(1, pos + downstream - 1)
        else:
            overhang = pos + downstream - length
            value = sub(pos - upstream, length) + sub(1, overhang - 1)
    if truncated and not cds and len(value) < 150:
        return ""
    return f">{tag}\n{value}\n"


@dataclass(frozen=True)
class TargetRegion:
    """A strand-oriented target with a 0-based, half-open prediction interval.

    ``sequence`` includes genomic accessibility context; ``requested`` excludes
    it. Instances are in-memory values and never modify their GenBank source.
    """

    tag: str
    sequence: str
    start: int
    end: int

    @property
    def requested(self) -> str:
        """Return only the requested target sequence, without accessibility flanks."""
        return self.sequence[self.start : self.end]


def _genomic_slice(sequence, start, end, circular):
    # Context may cross either end or more than one circular origin. Linear
    # records stop at their physical ends; never fabricate missing nucleotides.
    if circular:
        length = len(sequence)
        offset = start % length
        return (sequence * ((offset + end - start + length - 1) // length))[
            offset : offset + end - start
        ]
    return sequence[max(0, start) : min(len(sequence), end)]


def extract_targets(
    path: Path, up: int, down: int, region: str, context: int | None = 0
):
    """Yield requested CDS regions with additional genomic accessibility context.

    Preserve the established requested-region sequence, including spliced CDS
    order and omission of truncated UTRs shorter than 150 nt. Add ``context`` bases
    at both biological ends, wrapping circular records and clipping linear ones.
    None requests a full replicon length on each side for an unlimited span.
    TargetRegion coordinates are 0-based, half-open within the extended sequence.
    Stream compressed/plain GenBank without modifying it; reject negative context.
    """
    if context is not None and context < 0:
        raise ValueError("Accessibility context must be nonnegative")
    for record in genbank_records(path):
        source = str(record.seq)
        circular = record.annotations.get("topology") == "circular"
        for feature in record.features:
            if feature.type != "CDS":
                continue
            tag = feature.qualifiers.get("locus_tag", [""])[0]
            parts = list(feature.location.parts)
            # BioPerl sub_Location preserves GenBank location-list order.
            # Biopython reverses that list for a complemented compound location.
            if feature.location.strand == -1:
                parts.reverse()
            start, end = int(feature.location.start) + 1, int(feature.location.end)
            strand = 1 if feature.location.strand == 1 else -1
            split = len(parts) > 1

            def flank(pos, upstream, downstream, cds=False):
                return get_region(
                    source, circular, strand, pos, upstream, downstream, tag, cds=cds
                )

            if region == "5utr":
                pos = int(parts[0].start) + 1 if strand == 1 else int(parts[-1].end)
                requested = flank(pos, up, down)
                left, right = (
                    (pos - up - 1, pos + down - 1)
                    if strand == 1
                    else (pos + up, pos - down)
                )
            elif region == "3utr":
                pos = int(parts[-1].end) + 1 if strand == 1 else int(parts[0].start)
                requested = flank(pos, up, down)
                left, right = (
                    (pos - up - 1, pos + down - 1)
                    if strand == 1
                    else (pos + up, pos - down)
                )
            elif region == "cds":
                protein_region = str(feature.extract(record.seq))
                # Perl spliced_seq is returned in biological order, independent of
                # the order of joined locations in an atypical GenBank annotation.
                if split:
                    before = flank(
                        int(parts[0].start) + 1 if strand == 1 else int(parts[-1].end),
                        up,
                        0,
                        True,
                    ).removesuffix("\n")
                    after = flank(
                        int(parts[-1].end) + 1 if strand == 1 else int(parts[0].start),
                        0,
                        down,
                        True,
                    ).split("\n", 1)[1]
                elif strand == 1:
                    before = flank(
                        start + int(start == 1), up + int(start == 1), 0, True
                    ).removesuffix("\n")
                    if start == 1:
                        before = before[:-1]
                    after = flank(end, 0, down + 1, True).split("\n", 1)[1][1:]
                else:
                    before = flank(end - 1, up + 1, 0, True).removesuffix("\n")[:-1]
                    after = flank(
                        start if start == 1 else start - 1,
                        0,
                        down + int(start == 1),
                        True,
                    ).split("\n", 1)[1]
                    if start == 1:
                        after = after[1:]
                requested = before + protein_region + after
                left, right = (
                    (int(parts[0].start) - up, int(parts[-1].end) + down)
                    if strand == 1
                    else (int(parts[-1].end) + up, int(parts[0].start) - down)
                )
            else:
                raise ValueError(f"Unknown target region: {region}")
            if not requested:
                continue
            sequence = requested.split("\n", 1)[1].rstrip("\n")
            if not circular:
                left, right = (min(len(source), max(0, p)) for p in (left, right))
            amount = len(source) if context is None else context
            if strand == 1:
                before = _genomic_slice(source, left - amount, left, circular)
                after = _genomic_slice(source, right, right + amount, circular)
            else:
                before = _reverse_complement(
                    _genomic_slice(source, left, left + amount, circular)
                )
                after = _reverse_complement(
                    _genomic_slice(source, right - amount, right, circular)
                )
            yield TargetRegion(
                tag, before + sequence + after, len(before), len(before) + len(sequence)
            )


def parse_regions(path: Path, up: int, down: int, region: str) -> str:
    """Return requested CDS FASTA regions without accessibility context.

    Preserve the established 1-based flank conventions, strand, circular wrap,
    spliced CDS and short-truncated-UTR exclusion. Stream compressed/plain
    GenBank records without writing files; invalid coordinates raise ValueError.
    """
    return "".join(
        f">{r.tag}\n{r.requested}\n" for r in extract_targets(path, up, down, region)
    )


def check_gene_cds(path: Path, display_name: str | None = None) -> str:
    """Return a diagnostic when genes are annotated but no CDS features exist."""
    genes = cds = 0
    for record in genbank_records(path):
        genes += sum(f.type == "gene" for f in record.features)
        cds += sum(f.type == "CDS" for f in record.features)
    if cds == 0 and genes != 0:
        return f"CDS features are {cds} while gene features are {genes} in {display_name or path}\n"
    return ""


def extract_proteins(path: Path, organism: str | None = None) -> str:
    """Return translated CDS FASTA with organism:locus-tag identifiers.

    Use an explicit assembly code to group chromosomes and plasmids in DomClust,
    or the file's accession for a standalone helper. No organism catalog is used.
    Reject whitespace/colon in codes; stream the input without writing files.
    """
    code = organism or accession(path).replace(".", "_")
    if not code or re.search(r"\s|:", code):
        raise ValueError(
            "Protein organism code must be nonempty without whitespace or colon"
        )
    output = []
    for record in genbank_records(path):
        for feature in record.features:
            q = feature.qualifiers
            if feature.type == "CDS" and "locus_tag" in q and "translation" in q:
                output.append(f">{code}:{q['locus_tag'][0]}\n{q['translation'][0]}\n")
    return "".join(output)


def extract_16s(groups: list[str], directory: Path) -> str:
    """Return one annotated 16S sequence per organism, searching replicons in order."""
    output = []
    for group in groups:
        filenames = group.split(",")
        label = accession(filenames[0])
        found = False
        for filename in filenames:
            for record in genbank_records(directory / filename):
                for feature in record.features:
                    product = feature.qualifiers.get("product", [""])[0]
                    if feature.type == "rRNA" and re.search(
                        r"16S|Small subunit ribosomal RNA", product, re.I
                    ):
                        output.append(f">{label}\n{feature.extract(record.seq)}\n")
                        found = True
                        break
                if found:
                    break
            if found:
                break
    return "".join(output)
