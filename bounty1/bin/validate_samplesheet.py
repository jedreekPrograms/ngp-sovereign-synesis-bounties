#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

CONDITIONS = {'WT', *(f'SIRT{i}' for i in range(1, 8))}
MARKS = {'H3K9ac', 'H3K56ac'}


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
            matrix[(cond, rep)].add('WGBS')
        else:
            raise ValueError(f'unknown assay {assay}')
        if not r['fastq_1']:
            raise ValueError(f'{r["sample_id"]} is missing fastq_1')

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
