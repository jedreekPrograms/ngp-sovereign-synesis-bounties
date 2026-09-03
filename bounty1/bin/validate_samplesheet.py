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
REQUIRED_COLUMNS = {'sample_id','condition','replicate','assay','mark','fastq_1','fastq_2'}


def _read_rows(path: Path):
    with path.open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f'{path} is empty')
    missing = REQUIRED_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f'{path} missing columns: {sorted(missing)}')
    ids = [r['sample_id'] for r in rows]
    dup = [k for k,v in Counter(ids).items() if v > 1]
    if dup:
        raise ValueError(f'duplicate sample_id values: {dup}')
    return rows


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


def _validate_integrity_fields(row):
    _validate_fastq_url(row['sample_id'], row['fastq_1'], 'fastq_1')
    if row['fastq_2']:
        _validate_fastq_url(row['sample_id'], row['fastq_2'], 'fastq_2')
    _validate_optional_accessions(row)

    run = row.get('run_accession', '')
    if run:
        expected_prefix = f'https://download.cncb.ac.cn/gsa-human/HRA003336/{run}/'
        for fastq_field in ('fastq_1', 'fastq_2'):
            value = row.get(fastq_field, '')
            if value.startswith('https://download.cncb.ac.cn/') and not value.startswith(expected_prefix):
                raise ValueError(
                    f'{row["sample_id"]} {fastq_field} does not match run accession {run}: {value}'
                )
    return run


def validate(path: Path):
    rows = _read_rows(path)
    matrix = defaultdict(set)
    runs = []
    for r in rows:
        cond = r['condition']
        if cond not in CONDITIONS:
            raise ValueError(f'unknown condition {cond}')
        rep = int(r['replicate'])
        assay = r['assay']
        if assay == 'CHIP':
            if rep < 1:
                raise ValueError('CHIP replicate must be >= 1')
            mark = r['mark']
            if mark not in MARKS:
                raise ValueError(f'CHIP sample {r["sample_id"]} has invalid mark {mark}')
            matrix[(cond, rep)].add(mark)
        elif assay == 'WGBS':
            if rep < 1:
                raise ValueError('WGBS replicate must be >= 1')
            if r['mark']:
                raise ValueError(f'WGBS sample {r["sample_id"]} must have an empty mark')
            matrix[(cond, rep)].add('WGBS')
        else:
            raise ValueError(f'unknown assay {assay}; primary samplesheet only accepts CHIP/WGBS')

        run = _validate_integrity_fields(r)
        if run:
            runs.append(run)

    dup_runs = [k for k,v in Counter(runs).items() if v > 1]
    if dup_runs:
        raise ValueError(f'duplicate run accessions: {dup_runs}')

    paired = [key for key, assays in matrix.items() if {'H3K9ac','H3K56ac','WGBS'} <= assays]
    if len(paired) < 3:
        raise ValueError(f'need >=3 condition/replicate keys with H3K9ac, H3K56ac and WGBS; found {len(paired)}')
    return len(rows), len(paired)


def validate_controls(path: Path):
    rows = _read_rows(path)
    conditions = set()
    runs = []
    for r in rows:
        cond = r['condition']
        if cond not in CONDITIONS:
            raise ValueError(f'unknown INPUT condition {cond}')
        if r['assay'] != 'INPUT':
            raise ValueError(f'{r["sample_id"]} must have assay INPUT')
        if int(r['replicate']) != 0:
            raise ValueError(f'INPUT sample {r["sample_id"]} must use replicate 0')
        if r['mark']:
            raise ValueError(f'INPUT sample {r["sample_id"]} must have an empty mark')
        if cond in conditions:
            raise ValueError(f'duplicate INPUT control for {cond}')
        conditions.add(cond)
        run = _validate_integrity_fields(r)
        if run:
            runs.append(run)

    if conditions != CONDITIONS:
        missing = sorted(CONDITIONS - conditions)
        extra = sorted(conditions - CONDITIONS)
        raise ValueError(f'INPUT control conditions mismatch; missing={missing}, extra={extra}')
    if len(runs) != len(set(runs)):
        raise ValueError('duplicate run accession in INPUT controls')
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('samplesheet', type=Path)
    ap.add_argument('--controls-sheet', type=Path)
    args = ap.parse_args()
    n, paired = validate(args.samplesheet)
    msg = f'OK: {n} primary rows; {paired} fully paired condition/replicate keys'
    if args.controls_sheet:
        controls = validate_controls(args.controls_sheet)
        msg += f'; {controls} per-condition INPUT controls'
    print(msg)


if __name__ == '__main__':
    main()
