"""Validated deterministic statistical preprocessing from the legacy R scripts.

GEV fitting and correlated p-value joining have separate migration requirements;
no substitute numerical algorithm is used here.
"""
from collections import Counter
import math
from pathlib import Path
import re


def bh_adjust(pvalues: list[float]) -> list[float]:
    valid = [i for i,value in enumerate(pvalues) if not math.isnan(value)]
    ordered = sorted(valid,key=lambda i:pvalues[i],reverse=True)
    adjusted = [math.nan]*len(pvalues)
    current = 1.
    for rank,index in enumerate(ordered):
        current = min(current,len(valid)/(len(valid)-rank)*pvalues[index])
        adjusted[index] = current
    return adjusted


def r_number(value: float) -> str:
    if math.isnan(value):
        return 'NA'
    if math.isinf(value):
        return 'Inf' if value>0 else '-Inf'
    return format(value,'.15g')


def adjust_pvalues(path: Path, directory: Path) -> str:
    rows = [line.split() for line in path.read_text().splitlines()[1:] if line.strip()]
    width = len(rows[0]) if rows else 0
    values = [float(row[col]) if row[col] != 'NA' else math.nan for col in range(width) for row in rows]
    adjusted = bh_adjust(values)
    (directory/'padj.csv').write_text('fdr\n'+''.join(r_number(p)+'\n' for p in adjusted))
    return ''


def prepare_copra2(path: Path, directory: Path) -> str:
    rows = [line.split(';') for line in path.read_text().splitlines()[1:]]
    groups = dict.fromkeys(row[-1] for row in rows)
    output = []
    for cluster in groups:
        members = sorted((row for row in rows if row[-1]==cluster),key=lambda row:row[1])
        output.append('0;'+';'.join(row[0] for row in members)+';\n')
    (directory/'CopraRNA2_prep.csv').write_text(''.join(output))
    return ''


def filter_clusters(directory: Path) -> str:
    options = (directory/'CopraRNA_option_file.txt').read_text().splitlines()
    relative_size = float(options[4].split(':')[1])
    organism_count = sum('>' in line.split()[0] for line in (directory/'16s_sequences.fa').read_text().splitlines() if line.split())
    lines = (directory/'opt_tags.clustered').read_text().splitlines()
    header = lines[0].split(';')
    # read.table applies make.names; the wrong legacy rename is intentional.
    header = [name.replace('-','.') for name in header]
    header[34] = 'p-value'
    rows = [line.split(';') for line in lines[1:]]
    counts = Counter(int(row[-1]) for row in rows)
    last_cluster = int(rows[-1][-1])
    keep = {cluster for cluster,count in counts.items() if 1<=cluster<=last_cluster and count/organism_count>=relative_size}
    output = [';'.join(header)+'\n']
    output += [';'.join(row)+'\n' for row in rows if int(row[-1]) in keep]
    (directory/'opt_tags.clustered_rcsize').write_text(''.join(output))
    return ''
