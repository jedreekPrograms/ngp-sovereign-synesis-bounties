#!/usr/bin/env python3
"""Convert MethylDackel CpG bedGraph to Bismark coverage coordinates.

MethylDackel emits 0-based half-open bedGraph coordinates. Bismark coverage
files use 1-based CpG positions. The downstream DunedinPACE projector consumes
the latter convention, so this adapter performs only the coordinate conversion
and preserves the measured methylated/unmethylated counts.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def convert_line(line: str) -> str | None:
    """Convert one MethylDackel data row; ignore track/comment/blank rows."""
    stripped = line.strip()
    if not stripped or stripped.startswith(("#", "track ", "browser ")):
        return None
    fields = stripped.split("\t")
    if len(fields) < 6:
        raise ValueError(f"Expected at least 6 tab-separated columns: {stripped}")

    chrom = fields[0]
    start0 = int(fields[1])
    end0 = int(fields[2])
    percent = float(fields[3])
    methylated = int(fields[4])
    unmethylated = int(fields[5])
    if start0 < 0 or end0 <= start0:
        raise ValueError(f"Invalid bedGraph interval: {chrom}:{start0}-{end0}")
    if not 0.0 <= percent <= 100.0:
        raise ValueError(f"Invalid methylation percentage: {percent}")
    if methylated < 0 or unmethylated < 0:
        raise ValueError("Methylation counts must be non-negative")

    pos1 = start0 + 1
    return f"{chrom}\t{pos1}\t{pos1}\t{percent:g}\t{methylated}\t{unmethylated}\n"


def convert_file(source: Path, destination: Path) -> int:
    """Convert a complete MethylDackel output file to gzip Bismark coverage."""
    rows = 0
    with source.open("rt", encoding="utf-8") as src, gzip.open(
        destination, "wt", encoding="utf-8"
    ) as dst:
        for line in src:
            converted = convert_line(line)
            if converted is not None:
                dst.write(converted)
                rows += 1
    if rows == 0:
        raise ValueError(f"No methylation rows found in {source}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = convert_file(args.input, args.output)
    print(f"converted_rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
