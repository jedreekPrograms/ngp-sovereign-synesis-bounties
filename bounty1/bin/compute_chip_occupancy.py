#!/usr/bin/env python3
"""Compute SIRT6-target-locus and global H3K9ac/H3K56ac occupancy.

Primary bounty endpoint
-----------------------
The bounty asks for differential H3K9ac/H3K56ac occupancy at SIRT6 target
loci.  When --loci-bed is provided, the primary ``differential_occupancy``
column is therefore the fragment-level log2(CPM + 1) signal over an
independently defined SIRT6 CUT&RUN locus set, centered on the WT mean for the
same histone mark.  These loci are defined without using histone occupancy or
DunedinPACE values.

Secondary global endpoint
-------------------------
For each histone mark, a fixed peak universe is also created from the union of
all MACS3 peaks for that mark.  Each sample is quantified over the same
universe and reported as ``global_differential_occupancy``.  This is retained
as a prespecified sensitivity analysis and is never substituted for the SIRT6
locus endpoint after results are observed.

For paired-end libraries, only the primary read1 of each proper pair represents
a fragment. For single-end libraries, each primary mapped read is a fragment.
Query names are de-duplicated across regions so a fragment spanning more than
one interval is counted once.
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
    match = SAMPLE_RE.search(Path(name).name)
    if not match:
        raise ValueError(f'cannot parse sample metadata from {name}')
    return match.group(1), int(match.group(2)), match.group(3)


def merge_intervals(intervals):
    by_chr = defaultdict(list)
    for chrom, start, end in intervals:
        if end <= start:
            continue
        by_chr[chrom].append((start, end))
    merged = {}
    for chrom, values in by_chr.items():
        values.sort()
        out = []
        for start, end in values:
            if not out or start > out[-1][1]:
                out.append([start, end])
            else:
                out[-1][1] = max(out[-1][1], end)
        merged[chrom] = [(start, end) for start, end in out]
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
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f'no peak files matched {pattern}')
    for filename in files:
        _, _, mark = parse_sample(filename)
        with open(filename, encoding='utf-8') as handle:
            for line in handle:
                if line.strip():
                    chrom, start, end, *_ = line.rstrip().split('\t')
                    peaks_by_mark[mark].append((chrom, int(start), int(end)))
    missing = [mark for mark in MARKS if not peaks_by_mark[mark]]
    if missing:
        raise ValueError(f'missing peak files/intervals for: {", ".join(missing)}')
    return {mark: merge_intervals(peaks_by_mark[mark]) for mark in MARKS}


def is_fragment_representative(read):
    """Return True exactly once per usable DNA fragment."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.is_paired:
        return read.is_proper_pair and read.is_read1
    return True


def total_fragments(bam):
    return sum(1 for read in bam.fetch(until_eof=True) if is_fragment_representative(read))


def fragments_in_regions(bam, regions):
    names = set()
    for chrom, intervals in regions.items():
        if chrom not in bam.references:
            continue
        for start, end in intervals:
            for read in bam.fetch(chrom, start, end):
                if is_fragment_representative(read):
                    names.add(read.query_name)
    return len(names)


def log2_cpm(numerator, denominator):
    if denominator < 1:
        raise ValueError('fragment denominator must be positive')
    return math.log2((numerator / denominator * 1_000_000.0) + 1.0)


def center_on_wt(rows, value_key, output_key):
    means = {}
    for mark in MARKS:
        values = [
            row[value_key]
            for row in rows
            if row['condition'] == 'WT' and row['mark'] == mark and row[value_key] is not None
        ]
        if not values:
            raise ValueError(f'missing WT controls for {mark} / {value_key}')
        means[mark] = sum(values) / len(values)
    for row in rows:
        value = row[value_key]
        row[output_key] = None if value is None else value - means[row['mark']]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--peaks-glob', required=True)
    parser.add_argument('--bam-glob', required=True)
    parser.add_argument('--loci-bed', default='')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    try:
        import pysam
    except ImportError as exc:
        raise SystemExit('pysam is required for production occupancy calculation') from exc

    global_regions = peak_union_by_mark(args.peaks_glob)
    sirt6_regions = read_bed(args.loci_bed) if args.loci_bed else None

    rows = []
    bam_files = sorted(glob.glob(args.bam_glob))
    if not bam_files:
        raise ValueError(f'no BAM files matched {args.bam_glob}')

    for bam_filename in bam_files:
        condition, replicate, mark = parse_sample(bam_filename)
        with pysam.AlignmentFile(bam_filename, 'rb') as bam:
            denominator = total_fragments(bam)
        if denominator < 1:
            raise ValueError(f'no usable mapped fragments in {bam_filename}')

        with pysam.AlignmentFile(bam_filename, 'rb') as bam:
            global_count = fragments_in_regions(bam, global_regions[mark])
        global_value = log2_cpm(global_count, denominator)

        sirt6_value = None
        sirt6_count = None
        if sirt6_regions is not None:
            with pysam.AlignmentFile(bam_filename, 'rb') as bam:
                sirt6_count = fragments_in_regions(bam, sirt6_regions)
            sirt6_value = log2_cpm(sirt6_count, denominator)

        rows.append({
            'condition': condition,
            'replicate': replicate,
            'mark': mark,
            'total_fragments': denominator,
            'global_fragments': global_count,
            'global_log2_cpm': global_value,
            'sirt6_fragments': sirt6_count,
            'sirt6_log2_cpm': sirt6_value,
        })

    # Keep the global fixed-universe endpoint as a prespecified sensitivity
    # analysis, but make SIRT6 target loci the acceptance-facing endpoint.
    center_on_wt(rows, 'global_log2_cpm', 'global_differential_occupancy')
    if sirt6_regions is not None:
        center_on_wt(rows, 'sirt6_log2_cpm', 'sirt6_differential_occupancy')
        for row in rows:
            row['differential_occupancy'] = row['sirt6_differential_occupancy']
    else:
        for row in rows:
            row['sirt6_differential_occupancy'] = None
            row['differential_occupancy'] = row['global_differential_occupancy']

    fieldnames = [
        'condition', 'replicate', 'mark',
        'differential_occupancy',
        'global_differential_occupancy', 'global_log2_cpm', 'global_fragments',
        'sirt6_differential_occupancy', 'sirt6_log2_cpm', 'sirt6_fragments',
        'total_fragments',
    ]
    with open(args.output, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
