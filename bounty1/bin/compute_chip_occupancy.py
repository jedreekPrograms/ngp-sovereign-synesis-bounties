#!/usr/bin/env python3
"""Compute a pre-specified differential occupancy scalar per sample/mark.

Production execution uses pysam to count mapped fragments in the union of MACS3
peaks for each mark. The scalar is log2(CPM + 1) minus the WT mean for the same
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


def parse_sample(name):
    m = SAMPLE_RE.search(Path(name).name)
    if not m:
        raise ValueError(f'cannot parse sample metadata from {name}')
    return m.group(1), int(m.group(2)), m.group(3)


def merge_intervals(intervals):
    by_chr = defaultdict(list)
    for chrom, start, end in intervals:
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
        merged[chrom] = out
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--peaks-glob', required=True)
    ap.add_argument('--bam-glob', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    try:
        import pysam
    except ImportError as e:
        raise SystemExit('pysam is required for production occupancy calculation') from e

    peaks_by_mark = defaultdict(list)
    for fn in glob.glob(a.peaks_glob):
        _, _, mark = parse_sample(fn)
        with open(fn, encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    c, s, e, *_ = line.rstrip().split('\t')
                    peaks_by_mark[mark].append((c, int(s), int(e)))
    union = {m: merge_intervals(v) for m, v in peaks_by_mark.items()}

    raw = []
    for bam_fn in glob.glob(a.bam_glob):
        cond, rep, mark = parse_sample(bam_fn)
        bam = pysam.AlignmentFile(bam_fn, 'rb')
        count = 0
        for chrom, intervals in union[mark].items():
            if chrom not in bam.references:
                continue
            for s, e in intervals:
                count += bam.count(contig=chrom, start=s, end=e, read_callback='all')
        mapped = max(bam.mapped, 1)
        cpm = count / mapped * 1_000_000.0
        raw.append({'condition': cond, 'replicate': rep, 'mark': mark, 'log2_cpm': math.log2(cpm + 1.0)})
        bam.close()

    wt_mean = {}
    for mark in ('H3K9ac','H3K56ac'):
        vals = [r['log2_cpm'] for r in raw if r['condition']=='WT' and r['mark']==mark]
        if not vals:
            raise ValueError(f'missing WT controls for {mark}')
        wt_mean[mark] = sum(vals)/len(vals)
    for r in raw:
        r['differential_occupancy'] = r['log2_cpm'] - wt_mean[r['mark']]

    with open(a.output, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['condition','replicate','mark','differential_occupancy','log2_cpm'])
        w.writeheader(); w.writerows(raw)

if __name__ == '__main__':
    main()
