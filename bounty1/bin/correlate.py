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

MARKS = ('H3K9ac', 'H3K56ac')


def _pearson(x, y):
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('at least 3 paired observations are required')
    if pearsonr is not None:
        r, p = pearsonr(x, y)
        return float(r), float(p)
    mx, my = mean(x), mean(y)
    dx = [value - mx for value in x]
    dy = [value - my for value in y]
    denominator = math.sqrt(sum(value * value for value in dx) * sum(value * value for value in dy))
    if denominator == 0:
        raise ValueError('zero variance in one of the variables')
    r = sum(a * b for a, b in zip(dx, dy)) / denominator
    return r, float('nan')


def _key(row):
    return row['condition'], int(row['replicate'])


def read_pace(path):
    with open(path, newline='', encoding='utf-8') as handle:
        return {_key(row): float(row['dunedinpace']) for row in csv.DictReader(handle)}


def read_occ(path, column='differential_occupancy'):
    result = {}
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if column not in (reader.fieldnames or []):
            raise ValueError(f'occupancy column {column!r} not present in {path}')
        for row in reader:
            value = row.get(column, '')
            if value in ('', None):
                continue
            key = _key(row)
            result.setdefault(key, {})[row['mark']] = float(value)
    return result


def compute(pace, occ):
    paired = []
    for key in sorted(set(pace) & set(occ)):
        marks = occ[key]
        if all(mark in marks for mark in MARKS):
            paired.append((key, pace[key], marks['H3K9ac'], marks['H3K56ac']))
    if len(paired) < 3:
        raise ValueError(f'need >=3 fully paired samples, found {len(paired)}')

    pace_values = [item[1] for item in paired]
    output = {
        'n_paired': len(paired),
        'samples': [f'{key[0]}_rep{key[1]}' for key, *_ in paired],
    }
    for index, mark in ((2, 'H3K9ac'), (3, 'H3K56ac')):
        r_value, p_value = _pearson([item[index] for item in paired], pace_values)
        output[f'{mark}_vs_DunedinPACE'] = {
            'pearson_r': r_value,
            'p_value': p_value,
            'n': len(paired),
        }
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pace', required=True)
    parser.add_argument('--occupancy', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    pace = read_pace(args.pace)
    primary = compute(pace, read_occ(args.occupancy, 'differential_occupancy'))
    result = {
        **primary,
        'primary_endpoint': {
            'occupancy_column': 'differential_occupancy',
            'definition': 'global fixed histone peak-universe log2 fragment-CPM centered on WT mean',
        },
    }

    # Report the independent SIRT6-locus analysis as a sensitivity analysis.
    # It is never substituted for the primary endpoint after seeing results.
    try:
        secondary = compute(pace, read_occ(args.occupancy, 'sirt6_differential_occupancy'))
    except ValueError:
        secondary = None
    if secondary is not None:
        result['secondary_sirt6_locus_analysis'] = secondary

    Path(args.output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
