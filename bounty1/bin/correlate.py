#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean

try:
    from scipy.stats import pearsonr
except Exception:  # pragma: no cover
    pearsonr = None


def _pearson(x, y):
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('at least 3 paired observations are required')
    if pearsonr is not None:
        r, p = pearsonr(x, y)
        return float(r), float(p)
    mx, my = mean(x), mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    den = math.sqrt(sum(v*v for v in dx) * sum(v*v for v in dy))
    if den == 0:
        raise ValueError('zero variance in one of the variables')
    r = sum(a*b for a, b in zip(dx, dy)) / den
    return r, float('nan')


def _key(row):
    return row['condition'], int(row['replicate'])


def read_pace(path):
    with open(path, newline='', encoding='utf-8') as f:
        return {_key(r): float(r['dunedinpace']) for r in csv.DictReader(f)}


def read_occ(path):
    result = {}
    with open(path, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            key = _key(r)
            result.setdefault(key, {})[r['mark']] = float(r['differential_occupancy'])
    return result


def compute(pace, occ):
    paired = []
    for key in sorted(set(pace) & set(occ)):
        marks = occ[key]
        if 'H3K9ac' in marks and 'H3K56ac' in marks:
            paired.append((key, pace[key], marks['H3K9ac'], marks['H3K56ac']))
    if len(paired) < 3:
        raise ValueError(f'need >=3 fully paired samples, found {len(paired)}')
    y = [x[1] for x in paired]
    out = {'n_paired': len(paired), 'samples': [f'{k[0]}_rep{k[1]}' for k, *_ in paired]}
    for idx, mark in [(2, 'H3K9ac'), (3, 'H3K56ac')]:
        r, p = _pearson([x[idx] for x in paired], y)
        out[f'{mark}_vs_DunedinPACE'] = {'pearson_r': r, 'p_value': p}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pace', required=True)
    ap.add_argument('--occupancy', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    result = compute(read_pace(a.pace), read_occ(a.occupancy))
    Path(a.output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
