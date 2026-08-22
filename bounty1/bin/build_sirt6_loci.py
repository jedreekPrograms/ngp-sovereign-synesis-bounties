#!/usr/bin/env python3
"""Build reproducible SIRT6-specific loci from HRA005392 CUT&RUN peaks.

The locus definition is independent of H3K9ac/H3K56ac:
1. require reciprocal overlap between the two WT SIRT6 CUT&RUN replicates;
2. build the consensus interval from the overlapping portion;
3. identify reproducible intervals in the two SIRT6-KO SIRT6-antibody replicates;
4. remove WT consensus loci that overlap those KO reproducible intervals.

This gives a conservative set of experimentally supported SIRT6 binding loci for
histone occupancy quantification.
"""
import argparse
import glob
import re
from collections import defaultdict
from pathlib import Path

NAME_RE = re.compile(r'^(WT|SIRT6_KO)_rep([12])_SIRT6_CUTRUN_peaks\.narrowPeak$')


def read_peaks(path):
    intervals = []
    with Path(path).open(encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                raise ValueError(f'{path}:{line_no}: expected narrowPeak/BED columns')
            start, end = int(fields[1]), int(fields[2])
            if end > start:
                intervals.append((fields[0], start, end))
    return intervals


def reciprocal_overlap(a, b):
    if a[0] != b[0]:
        return 0.0
    overlap = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if overlap == 0:
        return 0.0
    return min(overlap / (a[2] - a[1]), overlap / (b[2] - b[1]))


def reproducible_intersections(rep1, rep2, min_reciprocal=0.50):
    by_chr = defaultdict(lambda: [[], []])
    for interval in rep1:
        by_chr[interval[0]][0].append(interval)
    for interval in rep2:
        by_chr[interval[0]][1].append(interval)

    out = []
    for chrom in sorted(by_chr):
        left, right = by_chr[chrom]
        left.sort(key=lambda x: x[1])
        right.sort(key=lambda x: x[1])
        for a in left:
            for b in right:
                if b[1] >= a[2]:
                    break
                if b[2] <= a[1]:
                    continue
                if reciprocal_overlap(a, b) >= min_reciprocal:
                    start = max(a[1], b[1])
                    end = min(a[2], b[2])
                    if end > start:
                        out.append((chrom, start, end))
    return merge_intervals(out)


def merge_intervals(intervals):
    by_chr = defaultdict(list)
    for chrom, start, end in intervals:
        by_chr[chrom].append((start, end))
    out = []
    for chrom in sorted(by_chr):
        for start, end in sorted(by_chr[chrom]):
            if not out or out[-1][0] != chrom or start > out[-1][2]:
                out.append([chrom, start, end])
            else:
                out[-1][2] = max(out[-1][2], end)
    return [tuple(x) for x in out]


def overlaps_any(interval, other_intervals):
    chrom, start, end = interval
    for other_chrom, other_start, other_end in other_intervals:
        if other_chrom < chrom:
            continue
        if other_chrom > chrom:
            break
        if other_start >= end:
            break
        if other_end > start:
            return True
    return False


def build_loci(peak_files, min_reciprocal=0.50):
    groups = {}
    for path in peak_files:
        name = Path(path).name
        match = NAME_RE.match(name)
        if not match:
            continue
        groups[(match.group(1), int(match.group(2)))] = read_peaks(path)

    required = {
        ('WT', 1), ('WT', 2),
        ('SIRT6_KO', 1), ('SIRT6_KO', 2),
    }
    missing = sorted(required - set(groups))
    if missing:
        raise ValueError(f'missing SIRT6 CUT&RUN peak sets: {missing}')

    wt = reproducible_intersections(
        groups[('WT', 1)], groups[('WT', 2)], min_reciprocal,
    )
    ko = reproducible_intersections(
        groups[('SIRT6_KO', 1)], groups[('SIRT6_KO', 2)], min_reciprocal,
    )
    ko = sorted(ko)
    specific = [interval for interval in wt if not overlaps_any(interval, ko)]
    return merge_intervals(specific)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--peaks-glob', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-reciprocal-overlap', type=float, default=0.50)
    args = parser.parse_args()

    if not (0 < args.min_reciprocal_overlap <= 1):
        raise ValueError('--min-reciprocal-overlap must be in (0, 1]')

    peak_files = sorted(glob.glob(args.peaks_glob))
    loci = build_loci(peak_files, args.min_reciprocal_overlap)
    if not loci:
        raise ValueError('no SIRT6-specific reproducible loci remained after KO filtering')

    with Path(args.output).open('w', encoding='utf-8') as handle:
        for chrom, start, end in loci:
            handle.write(f'{chrom}\t{start}\t{end}\n')

    print(f'SIRT6-specific reproducible loci: {len(loci)}')


if __name__ == '__main__':
    main()
