#!/usr/bin/env python3
"""Split an interleaved FASTQ stream into paired FIFO/file outputs.

This helper deliberately performs no record transformation. It validates basic
FASTQ framing and matching read names, then writes R1/R2 records separately so
mappers that require two filenames can still consume a fully streaming fastp
output without materialising cleaned whole-library FASTQs on disk.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def canonical_name(header: bytes) -> bytes:
    token = header.split(maxsplit=1)[0]
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    return token


def read_record(stream) -> list[bytes] | None:
    first = stream.readline()
    if not first:
        return None
    record = [first, stream.readline(), stream.readline(), stream.readline()]
    if any(not line for line in record):
        raise SystemExit("truncated interleaved FASTQ record")
    if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
        raise SystemExit("invalid FASTQ framing in interleaved stream")
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1-output", required=True)
    parser.add_argument("--r2-output", required=True)
    args = parser.parse_args()

    count = 0
    # Open in this order consistently with the mapper invocation (R1, then R2).
    with Path(args.r1_output).open("wb", buffering=0) as r1_out, Path(
        args.r2_output
    ).open("wb", buffering=0) as r2_out:
        while True:
            r1 = read_record(sys.stdin.buffer)
            if r1 is None:
                break
            r2 = read_record(sys.stdin.buffer)
            if r2 is None:
                raise SystemExit("interleaved FASTQ ended after R1 without R2")
            if canonical_name(r1[0]) != canonical_name(r2[0]):
                raise SystemExit(
                    "interleaved FASTQ pair names differ: "
                    f"{r1[0].decode(errors='replace').strip()} vs "
                    f"{r2[0].decode(errors='replace').strip()}"
                )
            r1_out.writelines(r1)
            r2_out.writelines(r2)
            count += 1

    print(f"Split interleaved FASTQ pairs: {count}", file=sys.stderr)
    if count == 0:
        raise SystemExit("interleaved FASTQ contained zero pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
