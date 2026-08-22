#!/usr/bin/env python3
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

CONDITIONS = {'WT', *(f'SIRT{i}' for i in range(1, 8))}
MARKS = {'H3K9ac', 'H3K56ac'}
MD5_RE = re.compile(r'^[0-9a-fA-F]{32}$')
RUN_RE = re.compile(r'^HRR\d+$')
EXPERIMENT_RE = re.compile(r'^HRX\d+$')
SAMPLE_RE = re.compile(r'^HRS\d+$')


def _validate_optional_accessions(row):
    checks = (
        ('run_accession', RUN_RE),
        ('experiment_accession', EXPERIMENT_RE),
        ('gsa_sample_accession', SAMPLE_RE),
    )
    for field, pattern in checks:
        value = row.get(field, '')
        if value and not pattern.fullmatch(value):
            raise ValueError(f'{row["sample_id"]} has invalid {field}: {value}')

    for field in ('read1_md5', 'read2_md5'):
        value = row.get(field, '')
        if value and not MD5_RE.fullmatch(value):
            raise ValueError(f'{row["sample_id"]} has invalid {field}: {value}')


def _validate_fastq_url(sample_id, value, expected_suffix):
    if not value:
        raise ValueError(f'{sample_id} is missing {expected_suffix}')
    parsed = urlparse(value)
    if parsed.scheme in {'http', 'https', 'ftp'} and not parsed.netloc:
        raise ValueError(f'{sample_id} has malformed URL: {value}')
    if not value.endswith('.fq.gz') and not value.endswith('.fastq.gz'):
        raise ValueError(f'{sample_id} FASTQ does not end in .fq.gz/.fastq.gz: {value}')


def validate(path: Path):
    with path.open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    required = {'sample_id','condition','replicate','assay','mark','fastq_1','fastq_2'}
    if not rows:
        raise ValueError('samplesheet is empty')
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f'missing columns: {sorted(missing)}')

    ids = [r['sample_id'] for r in rows]
    dup = [k for k,v in Counter(ids).items() if v > 1]
    if dup:
        raise ValueError(f'duplicate sample_id values: {dup}')

    matrix = defaultdict(set)
    runs = []
    for r in rows:
        cond = r['condition']
        if cond not in CONDITIONS:
            raise ValueError(f'unknown condition {cond}')
        rep = int(r['replicate'])
        if rep < 1:
            raise ValueError('replicate must be >= 1')
        assay = r['assay']
        if assay == 'CHIP':
            mark = r['mark']
            if mark not in MARKS:
                raise ValueError(f'CHIP sample {r["sample_id"]} has invalid mark {mark}')
            matrix[(cond, rep)].add(mark)
        elif assay == 'WGBS':
            if r['mark']:
                raise ValueError(f'WGBS sample {r["sample_id"]} must have an empty mark')
            matrix[(cond, rep)].add('WGBS')
        else:
            raise ValueError(f'unknown assay {assay}')

        _validate_fastq_url(r['sample_id'], r['fastq_1'], 'fastq_1')
        if r['fastq_2']:
            _validate_fastq_url(r['sample_id'], r['fastq_2'], 'fastq_2')
        _validate_optional_accessions(r)

        run = r.get('run_accession', '')
        if run:
            runs.append(run)
            expected_prefix = f'https://download.cncb.ac.cn/gsa-human/HRA003336/{run}/'
            for fastq_field in ('fastq_1', 'fastq_2'):
                value = r.get(fastq_field, '')
                if value.startswith('https://download.cncb.ac.cn/') and not value.startswith(expected_prefix):
                    raise ValueError(
                        f'{r["sample_id"]} {fastq_field} does not match run accession {run}: {value}'
                    )

    dup_runs = [k for k,v in Counter(runs).items() if v > 1]
    if dup_runs:
        raise ValueError(f'duplicate run accessions: {dup_runs}')

    paired = [key for key, assays in matrix.items() if {'H3K9ac','H3K56ac','WGBS'} <= assays]
    if len(paired) < 3:
        raise ValueError(f'need >=3 condition/replicate keys with H3K9ac, H3K56ac and WGBS; found {len(paired)}')
    return len(rows), len(paired)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('samplesheet', type=Path)
    args = ap.parse_args()
    n, paired = validate(args.samplesheet)
    print(f'OK: {n} rows; {paired} fully paired condition/replicate keys')


if __name__ == '__main__':
    main()
