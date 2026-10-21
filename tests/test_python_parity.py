"""Every migrated fixture must match independently captured legacy output."""
import json
from pathlib import Path
import shutil

import pytest

from coprarna import fasta, genome, tables, phylogeny, annotation, statistics
from tests.support.legacy import ROOT, FIXTURES, canonicalize


def run_port(spec, directory):
    script = Path(spec['script']).name
    args = spec['args']
    if script == 'format_fasta.pl':
        return fasta.format_fasta(directory / args[0])
    if script == 'fasta2genefile.pl':
        return fasta.fasta_to_genefile(directory / args[0])
    if script == 'parse_region_from_genome.pl':
        return genome.parse_regions(directory / args[0], int(args[1]), int(args[2]), args[3])
    if script == 'check_for_gene_CDS_features.pl':
        return genome.check_gene_cds(directory / args[0], args[0])
    if script == 'get_CDS_from_gbk.pl':
        return genome.extract_proteins(directory / args[0], (ROOT / 'coprarna_aux/kegg2refseqnew.csv').read_text(), args[0])
    if script == 'parse_16s_from_gbk.pl':
        return genome.extract_16s(args, directory)
    if script == 'blast2homfile.pl':
        return tables.blast_to_homology(directory / args[0])
    if script == 'disentangle_sorted_intarna_CSV.pl':
        return tables.disentangle(directory)
    if script == 'cluster_intarna_csv.pl':
        return tables.cluster_predictions(directory)
    if script == 'get_amount_sampled_values_and_add_to_table.pl':
        return tables.add_sampled_count(directory / args[0], bool(int(args[1])))
    if script == 'parse_homologs_from_domclust_table.pl':
        return tables.additional_homologs(directory / args[0], directory / args[1])
    if script == 'make_enrichment_table_from_py_output.pl':
        return tables.enrichment_report(directory / args[0])
    if script == 'transform_distmat.pl':
        return phylogeny.transform_distance(directory / args[0])
    if script == 'compute_weights.pl':
        return phylogeny.compute_weights(directory / args[0], directory / args[1])
    if script == 'annotate_raw_output.pl':
        return annotation.annotate_raw(directory / args[0], directory / args[1], args[2:], directory)
    if script == 'add_GI_genename_annotation_intarna.pl':
        return annotation.annotate_intarna(directory)
    if script == 'calc_padj.R':
        return statistics.adjust_pvalues(directory / args[0], directory)
    if script == 'prep_cop2.R':
        return statistics.prepare_copra2(directory / args[0], directory)
    if script == 'remove_clusters_under_percantage.R':
        return statistics.filter_clusters(directory)
    raise NotImplementedError(script)


MIGRATED = {'format_fasta.pl', 'fasta2genefile.pl', 'parse_region_from_genome.pl',
            'check_for_gene_CDS_features.pl', 'get_CDS_from_gbk.pl', 'parse_16s_from_gbk.pl',
            'blast2homfile.pl', 'disentangle_sorted_intarna_CSV.pl', 'cluster_intarna_csv.pl',
            'get_amount_sampled_values_and_add_to_table.pl', 'parse_homologs_from_domclust_table.pl',
            'make_enrichment_table_from_py_output.pl', 'transform_distmat.pl', 'compute_weights.pl', 'add_GI_genename_annotation_intarna.pl',
            'calc_padj.R', 'prep_cop2.R', 'remove_clusters_under_percantage.R', 'annotate_raw_output.pl'}
CASES = [p.parent for p in sorted(FIXTURES.glob('*/case.json'))
         if Path(json.loads(p.read_text())['script']).name in MIGRATED]


@pytest.mark.parametrize('case', CASES, ids=[p.name for p in CASES])
def test_python_matches_golden_bytes(case, tmp_path):
    spec = json.loads((case / 'case.json').read_text())
    shutil.copytree(case / 'input', tmp_path, dirs_exist_ok=True)
    stdout = run_port(spec, tmp_path)
    actual = {'stdout': stdout.encode()}
    actual.update({name: (tmp_path / name).read_bytes() for name in spec['outputs']})
    for name, value in actual.items():
        expected = (case / 'expected' / name).read_bytes()
        assert canonicalize(name, value, spec) == canonicalize(name, expected, spec)
