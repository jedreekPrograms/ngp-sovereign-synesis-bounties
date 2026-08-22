import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from validate_samplesheet import validate


def test_validate_minimum_complete_pairs(tmp_path):
    p = tmp_path / 'samples.csv'
    fields = ['sample_id','condition','replicate','assay','mark','fastq_1','fastq_2']
    rows = []
    for cond in ['WT','SIRT1','SIRT2']:
        for assay, mark in [('CHIP','H3K9ac'),('CHIP','H3K56ac'),('WGBS','')]:
            rows.append({
                'sample_id': f'{cond}_rep1_{mark or assay}', 'condition': cond,
                'replicate': '1', 'assay': assay, 'mark': mark,
                'fastq_1': '/x/r1.fq.gz', 'fastq_2': '/x/r2.fq.gz'
            })
    with p.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    n, paired = validate(p)
    assert n == 9
    assert paired == 3


def test_hra003336_production_samplesheet_is_complete():
    samplesheet = ROOT / 'resources' / 'samplesheet.csv'
    n, paired = validate(samplesheet)
    assert n == 48
    assert paired == 16

    with samplesheet.open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))

    assert len({r['run_accession'] for r in rows}) == 48
    assert len({r['experiment_accession'] for r in rows}) == 48
    assert len({r['gsa_sample_accession'] for r in rows}) == 48

    expected = {
        (condition, replicate, assay)
        for condition in ['WT', *(f'SIRT{i}' for i in range(1, 8))]
        for replicate in ('1', '2')
        for assay in ('H3K9ac', 'H3K56ac', 'WGBS')
    }
    observed = {
        (r['condition'], r['replicate'], r['mark'] or 'WGBS')
        for r in rows
    }
    assert observed == expected

    for row in rows:
        run = row['run_accession']
        assert row['fastq_1'].endswith(f'/{run}_f1.fq.gz')
        assert row['fastq_2'].endswith(f'/{run}_r2.fq.gz')
