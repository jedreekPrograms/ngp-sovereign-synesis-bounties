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
