#!/usr/bin/env python3
"""Compute differential H3K9ac/H3K56ac occupancy per sample.

Production execution uses pysam to count mapped fragments in a pre-specified
set of loci. For the bounty analysis those loci are derived independently from
SIRT6 CUT&RUN (HRA005392), avoiding a circular definition based on the histone
marks being correlated with DunedinPACE.

If --loci-bed is omitted, the script falls back to the union of MACS3 peaks for
each histone mark. The scalar is log2(CPM + 1) minus the WT mean for the same
mark. No target correlation value is used at any point.
"""
import argparse
import csv
import glob
import math
import re
from collections import defaultdict
from pathlib import Path

SAMPLE_RE = re.compile(r'^(WT|SIRT[1-7])_rep([12])_(H3K9ac|H3K56ac)')
MARKS = ('H3K9ac', 'H3K56ac')


def parse_sample(name):
    m = SAMPLE_RE.search(Path(name).name)
    if not m:
        raise ValueError(f'cannot parse sample metadata from {name}')
    return m.group(1), int(m.group(2)), m.group(3)


def merge_intervals(intervals):
    by_chr = defaultdict(list)
    for chrom, start, end in intervals:
        if end <= start:
            continue
        by_chr[chrom].append((start, end))
    merged = {}
    for chrom, xs in by_chr.items():
        xs.sort()
        out = []
        for s, e in xs:
            if not out or s > out[-1][1]:
                out.append([s, e])
            else:
                out[-1][1] = max(out[-1][1], e)
        merged[chrom] = [(s, e) for s, e in out]
    return merged


def read_bed(path):
    intervals = []
    with Path(path).open(encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                raise ValueError(f'{path}:{line_no}: expected at least 3 BED columns')
            intervals.append((fields[0], int(fields[1]), int(fields[2])))
    if not intervals:
        raise ValueError(f'no loci found in {path}')
    return merge_intervals(intervals)


def peak_union_by_mark(pattern):
    peaks_by_mark = defaultdict(list)
    for fn in glob.glob(pattern):
        _, _, mark = parse_sample(fn)
        with open(fn, encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    c, s, e, *_ = line.rstrip().split('\t')
                    peaks_by_mark[mark].append((c, int(s), int(e)))
    missing = [mark for mark in MARKS if not peaks_by_mark[mark]]
    if missing:
        raise ValueError(f'missing peak files/intervals for: {", ".join(missing)}')
    return {mark: merge_intervals(peaks_by_mark[mark]) for mark in MARKS}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--peaks-glob', default='')
    ap.add_argument('--bam-glob', required=True)
    ap.add_argument('--loci-bed', default='')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    try:
        import pysam
    except ImportError as e:
        raise SystemExit('pysam is required for production occupancy calculation') from e

    if a.loci_bed:
        shared_loci = read_bed(a.loci_bed)
        regions_by_mark = {mark: shared_loci for mark in MARKS}
    else:
        if not a.peaks_glob:
            raise ValueError('either --loci-bed or --peaks-glob is required')
        regions_by_mark = peak_union_by_mark(a.peaks_glob)

    raw = []
    bam_files = sorted(glob.glob(a.bam_glob))
    if not bam_files:
        raise ValueError(f'no BAM files matched {a.bam_glob}')

    for bam_fn in bam_files:
        cond, rep, mark = parse_sample(bam_fn)
        with pysam.AlignmentFile(bam_fn, 'rb') as bam:
            count = 0
            for chrom, intervals in regions_by_mark[mark].items():
                if chrom not in bam.references:
                    continue
                for start, end in intervals:
                    count += bam.count(
                        contig=chrom,
                        start=start,
                        end=end,
                        read_callback='all',
                    )
            mapped = max(bam.mapped, 1)
            cpm = count / mapped * 1_000_000.0
            raw.append({
                'condition': cond,
                'replicate': rep,
                'mark': mark,
                'log2_cpm': math.log2(cpm + 1.0),
            })

    wt_mean = {}
    for mark in MARKS:
        vals = [
            row['log2_cpm'] for row in raw
            if row['condition'] == 'WT' and row['mark'] == mark
        ]
        if not vals:
            raise ValueError(f'missing WT controls for {mark}')
        wt_mean[mark] = sum(vals) / len(vals)

    for row in raw:
        row['differential_occupancy'] = row['log2_cpm'] - wt_mean[row['mark']]

    with open(a.output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                'condition', 'replicate', 'mark',
                'differential_occupancy', 'log2_cpm',
            ],
        )
        writer.writeheader()
        writer.writerows(raw)


if __name__ == '__main__':
    main()
