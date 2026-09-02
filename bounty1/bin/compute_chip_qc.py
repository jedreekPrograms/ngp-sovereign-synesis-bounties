#!/usr/bin/env python3
"""Compute measured ChIP-seq IDR and NRF quality metrics without fabrication.

NRF is calculated from coordinate-sorted, pre-deduplication ``*.position.bam``
checkpoints.  Shards for the same sample are merged before ``samtools markdup``
so duplicates crossing shard boundaries are counted correctly.  The reported
sample NRF is:

    non-duplicate primary mapped reads at MAPQ >= threshold
    -------------------------------------------------------
       all primary mapped reads at MAPQ >= threshold

IDR is calculated by the Kundaje lab ``idr`` implementation from the two
biological-replicate MACS3 narrowPeak files for each supplied condition/mark.
The mark-level IDR exposed to the organizer manifest is the worst measured
global IDR among peaks retained at the configured threshold, across all
supplied complete replicate pairs.  Pair-level counts and values are retained
in the JSON so the scalar acceptance field is auditable rather than a constant.

The script fails closed on missing replicate pairs, empty inputs, zero usable
reads, malformed IDR output, or a pair with zero reproducible peaks.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import re
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

MARKS = ("H3K9ac", "H3K56ac")
SAMPLE_RE = re.compile(r"(WT|SIRT[1-7])_rep([12])_(H3K9ac|H3K56ac)")
PRIMARY_EXCLUDE = 4 + 256 + 2048
PRIMARY_NON_DUP_EXCLUDE = PRIMARY_EXCLUDE + 1024


def sample_meta(path: str | Path) -> tuple[str, int, str, str]:
    name = Path(path).name
    match = SAMPLE_RE.search(name)
    if not match:
        raise ValueError(f"cannot parse ChIP sample metadata from {name}")
    condition, replicate, mark = match.group(1), int(match.group(2)), match.group(3)
    return condition, replicate, mark, f"{condition}_rep{replicate}_{mark}"


def run_count(command: list[str]) -> int:
    text = subprocess.check_output(command, text=True).strip()
    value = int(text)
    if value < 0:
        raise ValueError(f"negative count from command: {' '.join(command)}")
    return value


def compute_sample_nrf(
    sample_id: str,
    shards: list[Path],
    mapq: int,
    samtools: str,
    threads: int,
    workdir: Path,
) -> dict:
    if not shards:
        raise ValueError(f"no position BAM shards for {sample_id}")
    for shard in shards:
        if not shard.is_file() or shard.stat().st_size == 0:
            raise ValueError(f"missing/empty position BAM shard: {shard}")
        subprocess.run([samtools, "quickcheck", "-v", str(shard)], check=True)

    merged = workdir / f"{sample_id}.position.bam"
    marked = workdir / f"{sample_id}.marked.bam"
    if len(shards) == 1:
        subprocess.run([samtools, "view", "-b", "-o", str(merged), str(shards[0])], check=True)
    else:
        subprocess.run(
            [samtools, "merge", "-f", "-@", str(threads), str(merged), *map(str, shards)],
            check=True,
        )
    subprocess.run([samtools, "quickcheck", "-v", str(merged)], check=True)
    subprocess.run(
        [samtools, "markdup", "-@", str(threads), str(merged), str(marked)],
        check=True,
    )
    subprocess.run([samtools, "quickcheck", "-v", str(marked)], check=True)

    total = run_count(
        [samtools, "view", "-c", "-q", str(mapq), "-F", str(PRIMARY_EXCLUDE), str(marked)]
    )
    distinct = run_count(
        [
            samtools,
            "view",
            "-c",
            "-q",
            str(mapq),
            "-F",
            str(PRIMARY_NON_DUP_EXCLUDE),
            str(marked),
        ]
    )
    if total <= 0:
        raise ValueError(f"{sample_id}: zero primary mapped reads at MAPQ >= {mapq}")
    if distinct > total:
        raise ValueError(f"{sample_id}: non-duplicate count {distinct} exceeds total {total}")
    nrf = distinct / total
    if not (0.0 <= nrf <= 1.0):
        raise ValueError(f"{sample_id}: invalid NRF {nrf}")
    return {
        "sample_id": sample_id,
        "mapq_threshold": mapq,
        "total_primary_mapped_reads": total,
        "nonduplicate_primary_mapped_reads": distinct,
        "nrf": nrf,
        "position_bam_shards": [str(path) for path in shards],
    }


def count_peaks(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def parse_idr_output(path: Path, threshold: float) -> dict:
    global_idrs: list[float] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{line_no}: expected >=12 IDR columns")
            transformed = float(fields[11])
            if not math.isfinite(transformed):
                raise ValueError(f"{path}:{line_no}: non-finite global IDR score")
            global_idr = 10.0 ** (-transformed)
            if not (0.0 <= global_idr <= 1.0):
                raise ValueError(f"{path}:{line_no}: invalid global IDR {global_idr}")
            global_idrs.append(global_idr)
    if not global_idrs:
        raise ValueError(f"IDR produced no comparable peaks: {path}")
    passing = [value for value in global_idrs if value <= threshold]
    if not passing:
        raise ValueError(f"no peaks pass global IDR <= {threshold} in {path}")
    return {
        "compared_peak_count": len(global_idrs),
        "reproducible_peak_count": len(passing),
        "worst_retained_global_idr": max(passing),
        "best_global_idr": min(global_idrs),
    }


def run_pair_idr(
    condition: str,
    mark: str,
    rep1: Path,
    rep2: Path,
    threshold: float,
    idr_command: list[str],
    workdir: Path,
) -> dict:
    if count_peaks(rep1) == 0 or count_peaks(rep2) == 0:
        raise ValueError(f"{condition} {mark}: empty narrowPeak replicate")
    output = workdir / f"{condition}_{mark}.idr.narrowPeak"
    command = [
        *idr_command,
        "--samples",
        str(rep1),
        str(rep2),
        "--input-file-type",
        "narrowPeak",
        "--rank",
        "p.value",
        "--soft-idr-threshold",
        str(threshold),
        "--output-file",
        str(output),
    ]
    subprocess.run(command, check=True)
    if not output.is_file() or output.stat().st_size == 0:
        raise ValueError(f"{condition} {mark}: IDR output missing/empty")
    summary = parse_idr_output(output, threshold)
    denominator = min(count_peaks(rep1), count_peaks(rep2))
    summary.update(
        {
            "condition": condition,
            "mark": mark,
            "replicate_1_peak_file": str(rep1),
            "replicate_2_peak_file": str(rep2),
            "minimum_input_peak_count": denominator,
            "reproducible_fraction_of_smaller_peak_set": (
                summary["reproducible_peak_count"] / denominator if denominator else 0.0
            ),
        }
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--position-bam-glob", required=True)
    parser.add_argument("--peaks-glob", required=True)
    parser.add_argument("--mapq", type=int, default=30)
    parser.add_argument("--idr-threshold", type=float, default=0.05)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--idr-command", default="idr")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if args.mapq < 0:
        raise ValueError("MAPQ must be non-negative")
    if not (0.0 < args.idr_threshold < 1.0):
        raise ValueError("IDR threshold must be between 0 and 1")
    if args.threads < 1:
        raise ValueError("threads must be >=1")

    bam_paths = [Path(p) for p in sorted(glob.glob(args.position_bam_glob, recursive=True))]
    peak_paths = [Path(p) for p in sorted(glob.glob(args.peaks_glob, recursive=True))]
    if not bam_paths:
        raise ValueError(f"no position BAMs matched {args.position_bam_glob}")
    if not peak_paths:
        raise ValueError(f"no narrowPeak files matched {args.peaks_glob}")

    bam_groups: dict[str, list[Path]] = defaultdict(list)
    sample_dims: dict[str, tuple[str, int, str]] = {}
    for path in bam_paths:
        condition, replicate, mark, sample_id = sample_meta(path)
        bam_groups[sample_id].append(path)
        sample_dims[sample_id] = (condition, replicate, mark)

    peak_groups: dict[tuple[str, str], dict[int, Path]] = defaultdict(dict)
    for path in peak_paths:
        condition, replicate, mark, _ = sample_meta(path)
        key = (condition, mark)
        if replicate in peak_groups[key]:
            raise ValueError(f"duplicate peak file for {condition} rep{replicate} {mark}")
        peak_groups[key][replicate] = path

    idr_command = shlex.split(args.idr_command)
    if not idr_command:
        raise ValueError("empty IDR command")

    with tempfile.TemporaryDirectory(prefix="bounty1-chip-qc-") as temp:
        workdir = Path(temp)
        sample_qc = []
        for sample_id, shards in sorted(bam_groups.items()):
            row = compute_sample_nrf(
                sample_id,
                sorted(shards),
                args.mapq,
                args.samtools,
                args.threads,
                workdir,
            )
            condition, replicate, mark = sample_dims[sample_id]
            row.update({"condition": condition, "replicate": replicate, "mark": mark})
            sample_qc.append(row)

        pair_qc = []
        for (condition, mark), reps in sorted(peak_groups.items()):
            if set(reps) != {1, 2}:
                raise ValueError(
                    f"{condition} {mark}: need exactly biological replicates 1 and 2, found {sorted(reps)}"
                )
            pair_qc.append(
                run_pair_idr(
                    condition,
                    mark,
                    reps[1],
                    reps[2],
                    args.idr_threshold,
                    idr_command,
                    workdir,
                )
            )

    result = {
        "schema_version": 1,
        "mapq_threshold": args.mapq,
        "idr_threshold": args.idr_threshold,
        "nrf_definition": (
            "non-duplicate primary mapped reads / all primary mapped reads at the configured MAPQ, "
            "after merging all pre-deduplication position-BAM shards for a sample"
        ),
        "idr_definition": (
            "Kundaje IDR 2.0.4.2 on biological-replicate MACS3 narrowPeak files ranked by p.value; "
            "mark-level idr is the worst measured global IDR among threshold-retained peaks across "
            "all supplied complete condition/replicate pairs"
        ),
    }

    for mark in MARKS:
        mark_samples = [row for row in sample_qc if row["mark"] == mark]
        mark_pairs = [row for row in pair_qc if row["mark"] == mark]
        if not mark_samples:
            raise ValueError(f"no pre-deduplication BAM evidence supplied for {mark}")
        if not mark_pairs:
            raise ValueError(f"no complete IDR replicate pair supplied for {mark}")
        result[mark] = {
            "idr": max(row["worst_retained_global_idr"] for row in mark_pairs),
            "nrf": min(row["nrf"] for row in mark_samples),
            "samples": mark_samples,
            "replicate_pairs": mark_pairs,
            "reproducible_peak_count": sum(row["reproducible_peak_count"] for row in mark_pairs),
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
