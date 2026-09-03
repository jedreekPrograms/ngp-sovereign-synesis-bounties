#!/usr/bin/env python3
"""Probe public archive FASTQ URLs and plan the smallest real-data pilot."""
import argparse
import csv
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def load_rows(paths):
    rows = []
    for path in paths:
        with Path(path).open(newline='', encoding='utf-8') as handle:
            for row in csv.DictReader(handle):
                row = dict(row)
                row['_sheet'] = str(path)
                rows.append(row)
    return rows


def probe_one(item, timeout, retries):
    row, mate = item
    url = row[mate]
    proc = subprocess.run(
        [
            'curl', '--head', '--location', '--silent', '--show-error',
            '--max-time', str(timeout), '--retry', str(retries), url,
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    codes = re.findall(r'(?im)^HTTP/\S+\s+(\d{3})', proc.stdout)
    lengths = [
        int(value)
        for value in re.findall(r'(?im)^content-length:\s*(\d+)\s*$', proc.stdout)
    ]
    return {
        'sample_id': row.get('sample_id', ''),
        'condition': row.get('condition', ''),
        'replicate': row.get('replicate', ''),
        'assay': row.get('assay', ''),
        'mark': row.get('mark') or row.get('target') or '',
        'mate': mate,
        'url': url,
        'status': int(codes[-1]) if codes else None,
        'bytes': lengths[-1] if lengths else None,
        'stderr': proc.stderr.strip(),
    }


def gib(value):
    return value / 1024**3


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sheet', action='append', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--timeout', type=int, default=25)
    parser.add_argument('--retries', type=int, default=2)
    args = parser.parse_args()

    rows = load_rows(args.sheet)
    work = []
    for row in rows:
        for mate in ('fastq_1', 'fastq_2'):
            if row.get(mate):
                work.append((row, mate))

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(probe_one, item, args.timeout, args.retries) for item in work]
        for future in as_completed(futures):
            results.append(future.result())

    results.sort(key=lambda row: (row['sample_id'], row['mate']))
    failures = [
        row for row in results
        if row['status'] is None or not 200 <= row['status'] < 400 or row['bytes'] is None
    ]

    with Path(args.output).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                'sample_id', 'condition', 'replicate', 'assay', 'mark',
                'mate', 'status', 'bytes', 'gib', 'url',
            ],
            delimiter='\t',
        )
        writer.writeheader()
        for row in results:
            writer.writerow({
                **{key: row[key] for key in writer.fieldnames if key not in {'gib'}},
                'gib': '' if row['bytes'] is None else f"{gib(row['bytes']):.3f}",
            })

    total = sum(row['bytes'] or 0 for row in results)
    known = sum(row['bytes'] is not None for row in results)
    print(f'FASTQ URLs checked: {len(results)}')
    print(f'HTTP/size failures: {len(failures)}')
    print(f'URLs exposing Content-Length: {known}/{len(results)}')
    print(f'Known compressed total: {gib(total):.2f} GiB')

    per_sample = defaultdict(int)
    for row in results:
        if row['bytes'] is not None:
            per_sample[row['sample_id']] += row['bytes']
    print('\nSmallest individual libraries:')
    for sample, size in sorted(per_sample.items(), key=lambda item: item[1])[:12]:
        print(f'  {sample:42s} {gib(size):8.2f} GiB')

    # A correlation observation needs H3K9ac, H3K56ac and WGBS from the same
    # condition + biological replicate. Rank complete observations by transfer
    # size so a pilot can use the smallest real paired data first.
    primary = defaultdict(lambda: {'marks': set(), 'bytes': 0})
    for row in results:
        if row['assay'] not in {'CHIP', 'WGBS'} or row['bytes'] is None:
            continue
        key = (row['condition'], row['replicate'])
        primary[key]['bytes'] += row['bytes']
        if row['assay'] == 'WGBS':
            primary[key]['marks'].add('WGBS')
        else:
            primary[key]['marks'].add(row['mark'])

    required = {'H3K9ac', 'H3K56ac', 'WGBS'}
    complete = [
        (key, value['bytes'])
        for key, value in primary.items()
        if required <= value['marks']
    ]
    complete.sort(key=lambda item: item[1])
    print('\nSmallest complete correlation observations (H3K9ac + H3K56ac + WGBS):')
    for (condition, replicate), size in complete[:10]:
        print(f'  {condition:8s} rep{replicate}: {gib(size):8.2f} GiB')
    if len(complete) >= 3:
        first_three = complete[:3]
        pilot_bytes = sum(size for _, size in first_three)
        labels = ', '.join(f'{cond}/rep{rep}' for (cond, rep), _ in first_three)
        print(f'\nSmallest 3-observation primary-data pilot: {labels} = {gib(pilot_bytes):.2f} GiB')

    if failures:
        for row in failures:
            print('FAIL:', row['sample_id'], row['mate'], row['status'], row['stderr'], row['url'])
        raise SystemExit(1)


if __name__ == '__main__':
    main()
