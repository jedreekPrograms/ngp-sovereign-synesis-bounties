#!/usr/bin/env python3
"""Pair-aware MAPQ filtering for name-sorted Bismark BAM files.

A paired fragment is retained only when exactly one primary R1 and one primary
R2 are present, both are mapped, and both satisfy the MAPQ threshold. This
prevents per-read MAPQ filtering from creating orphaned mates that break
paired-end Bismark methylation extraction.
"""
import argparse
import json
from itertools import groupby
from pathlib import Path


def usable_primary(read, mapq):
    return (
        not read.is_secondary
        and not read.is_supplementary
        and not read.is_unmapped
        and read.mapping_quality >= mapq
    )


def keep_paired_group(reads, mapq):
    primary = [read for read in reads if not read.is_secondary and not read.is_supplementary]
    if len(primary) != 2:
        return False
    read1 = [read for read in primary if read.is_read1]
    read2 = [read for read in primary if read.is_read2]
    if len(read1) != 1 or len(read2) != 1:
        return False
    return usable_primary(read1[0], mapq) and usable_primary(read2[0], mapq)


def filter_bam(input_path, output_path, mapq, paired=True, stats_path=''):
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit('pysam is required for BAM filtering') from exc

    stats = {
        'mapq_threshold': mapq,
        'paired': paired,
        'groups_seen': 0,
        'groups_kept': 0,
        'alignments_written': 0,
    }

    with pysam.AlignmentFile(input_path, 'rb') as source:
        header = source.header.to_dict()
        header.setdefault('HD', {})['SO'] = 'queryname'
        with pysam.AlignmentFile(output_path, 'wb', header=header) as target:
            if paired:
                for _, iterator in groupby(source.fetch(until_eof=True), key=lambda read: read.query_name):
                    reads = list(iterator)
                    stats['groups_seen'] += 1
                    if keep_paired_group(reads, mapq):
                        for read in reads:
                            if not read.is_secondary and not read.is_supplementary:
                                target.write(read)
                                stats['alignments_written'] += 1
                        stats['groups_kept'] += 1
            else:
                for read in source.fetch(until_eof=True):
                    stats['groups_seen'] += 1
                    if usable_primary(read, mapq):
                        target.write(read)
                        stats['groups_kept'] += 1
                        stats['alignments_written'] += 1

    if stats['alignments_written'] == 0:
        raise ValueError(f'no alignments passed pair-aware MAPQ >= {mapq}')
    if paired and stats['alignments_written'] != stats['groups_kept'] * 2:
        raise AssertionError('paired output did not contain exactly two primary alignments per kept group')

    if stats_path:
        Path(stats_path).write_text(json.dumps(stats, indent=2) + '\n', encoding='utf-8')
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--mapq', type=int, required=True)
    parser.add_argument('--paired', action='store_true')
    parser.add_argument('--stats', default='')
    args = parser.parse_args()
    stats = filter_bam(args.input, args.output, args.mapq, paired=args.paired, stats_path=args.stats)
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
