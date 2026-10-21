"""GenBank helpers preserving the legacy coordinate and boundary conventions."""
from pathlib import Path
import re

from Bio import SeqIO

# Perl's legacy transliteration complements only these four bases, not all IUPAC.
_COMPLEMENT = str.maketrans('gatcGATC', 'ctagCTAG')


def _reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(_COMPLEMENT)


def _subseq(sequence: str, start: int, end: int) -> str:
    # Bio::Seq allows an empty interval when start == end+1.
    if start < 1 or end > len(sequence) or start > end+1:
        raise ValueError(f'Invalid legacy subsequence coordinates: {start}, {end}')
    return sequence[start-1:end]


def get_region(sequence: str, circular: bool, strand: int, pos: int,
               upstream: int, downstream: int, tag: str, *, cds: bool = False) -> str:
    """Direct port of get_region/get_region_cds; positions are 1-based."""
    length = len(sequence)
    sub = lambda start, end: _subseq(sequence, start, end)
    truncated = False
    if strand == -1:
        if pos-downstream >= 1 and pos+upstream <= length:
            value = _reverse_complement(sub(pos-downstream+1, pos+upstream))
        elif not circular:
            truncated = True
            if pos-downstream < 1:
                value = _reverse_complement(sub(1, pos+upstream))
            else:
                value = _reverse_complement(sub(pos-downstream+1, length))
        elif pos-downstream < 1:
            overhang = downstream-pos
            value = _reverse_complement(sub(1, pos+upstream)) + _reverse_complement(sub(length-overhang+1, length))
        else:
            overhang = pos+upstream-length
            value = _reverse_complement(sub(1, overhang)) + _reverse_complement(sub(pos-downstream+1, length))
    else:
        if pos-upstream >= 1 and pos+downstream <= length:
            value = sub(pos-upstream, pos+downstream-1)
        elif not circular:
            truncated = True
            if pos-upstream < 1:
                value = sub(1, pos+downstream-1)
            else:
                value = sub(pos-upstream, length)
        elif pos-upstream < 1:
            overhang = upstream-pos
            value = sub(length-overhang, length) + sub(1, pos+downstream-1)
        else:
            overhang = pos+downstream-length
            value = sub(pos-upstream, length) + sub(1, overhang-1)
    if truncated and not cds and len(value) < 150:
        return ''
    return f'>{tag}\n{value}\n'


def parse_regions(path: Path, up: int, down: int, region: str) -> str:
    records = list(SeqIO.parse(path, 'genbank'))
    # Intentionally use the first sequence for features of every record, as Perl does.
    source = str(records[0].seq)
    circular = records[0].annotations.get('topology') == 'circular'
    output = []
    for record in records:
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            tag = feature.qualifiers.get('locus_tag', [''])[0]
            parts = list(feature.location.parts)
            # BioPerl sub_Location preserves GenBank location-list order.
            # Biopython reverses that list for a complemented compound location.
            if feature.location.strand == -1:
                parts.reverse()
            start, end = int(feature.location.start)+1, int(feature.location.end)
            strand = feature.location.strand
            split = len(parts) > 1
            def flank(pos, upstream, downstream, cds=False):
                return get_region(source, circular, strand, pos, upstream, downstream, tag, cds=cds)
            if region == '5utr':
                pos = int(parts[0].start)+1 if strand == 1 else int(parts[-1].end)
                output.append(flank(pos, up, down))
            elif region == '3utr':
                pos = int(parts[-1].end)+1 if strand == 1 else int(parts[0].start)
                output.append(flank(pos, up, down))
            elif region == 'cds':
                protein_region = str(feature.extract(record.seq))
                # Perl spliced_seq is returned in biological order, independent of
                # the order of joined locations in an atypical GenBank annotation.
                if split:
                    before = flank(int(parts[0].start)+1 if strand==1 else int(parts[-1].end), up, 0, True).removesuffix('\n')
                    after = flank(int(parts[-1].end)+1 if strand==1 else int(parts[0].start), 0, down, True).split('\n',1)[1]
                elif strand == 1:
                    before = flank(start+int(start==1), up+int(start==1), 0, True).removesuffix('\n')
                    if start == 1:
                        before = before[:-1]
                    after = flank(end, 0, down+1, True).split('\n',1)[1][1:]
                else:
                    before = flank(end-1, up+1, 0, True).removesuffix('\n')[:-1]
                    after = flank(start if start==1 else start-1, 0, down+int(start==1), True).split('\n',1)[1]
                    if start == 1:
                        after = after[1:]
                output.append(before + protein_region + after)
    return ''.join(output)


def check_gene_cds(path: Path, display_name: str | None = None) -> str:
    genes = cds = 0
    for record in SeqIO.parse(path, 'genbank'):
        genes += sum(f.type == 'gene' for f in record.features)
        cds += sum(f.type == 'CDS' for f in record.features)
    if cds == 0 and genes != 0:
        return f'CDS features are {cds} while gene features are {genes} in {display_name or path}\n'
    return ''


def extract_proteins(path: Path, mapping_text: str, display_name: str | None = None) -> str:
    # Matching is deliberately unanchored, and removes exactly three filename bytes.
    accession = (display_name or str(path))[:-3]
    code = ''
    for line in mapping_text.splitlines():
        if re.search(accession, line):
            code = re.split(r'\s', line)[0]
            break
    output = []
    for record in SeqIO.parse(path, 'genbank'):
        for feature in record.features:
            q = feature.qualifiers
            if feature.type == 'CDS' and 'locus_tag' in q and 'translation' in q:
                output.append(f'>{code}:{q["locus_tag"][0]}\n{q["translation"][0]}\n')
    return ''.join(output)


def extract_16s(groups: list[str], directory: Path) -> str:
    output = []
    for group in groups:
        filenames = group.split(',')
        label = filenames[0][:-3]
        found = False
        for filename in filenames:
            for record in SeqIO.parse(directory / filename, 'genbank'):
                for feature in record.features:
                    product = feature.qualifiers.get('product', [''])[0]
                    if feature.type == 'rRNA' and re.search(r'16S|Small subunit ribosomal RNA', product, re.I):
                        output.append(f'>{label}\n{feature.extract(record.seq)}\n')
                        found = True
                        break
                if found:
                    break
            if found:
                break
    return ''.join(output)
