#!/usr/bin/env python3
"""Bounded-disk Abismal screening for an interleaved cleaned FASTQ stream.

Abismal requires separate mate filenames. Feeding both mates through FIFOs can
 deadlock because the mapper is free to open/read its inputs sequentially while
an interleaved producer must advance R1 and R2 together. This helper avoids
that class of deadlock without materialising whole cleaned libraries: it writes
one bounded pair of temporary FASTQ chunks, maps that chunk, emits complete
candidate pairs to gzip outputs, then deletes the chunk before continuing.
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def canonical_name(header: bytes) -> bytes:
    token = header.split(maxsplit=1)[0]
    if token.startswith(b"@"):
        token = token[1:]
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    return token


def canonical_sam_name(name: str) -> str:
    return re.sub(r"/[12]$", "", name)


def read_record(stream) -> list[bytes] | None:
    first = stream.readline()
    if not first:
        return None
    record = [first, stream.readline(), stream.readline(), stream.readline()]
    if any(not line for line in record):
        raise RuntimeError("truncated interleaved FASTQ record")
    if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
        raise RuntimeError("invalid FASTQ framing in interleaved stream")
    return record


def write_chunk(stream, r1_path: Path, r2_path: Path, max_pairs: int) -> int:
    count = 0
    with r1_path.open("wb") as r1_out, r2_path.open("wb") as r2_out:
        while count < max_pairs:
            r1 = read_record(stream)
            if r1 is None:
                break
            r2 = read_record(stream)
            if r2 is None:
                raise RuntimeError("interleaved FASTQ ended after R1 without R2")
            if canonical_name(r1[0]) != canonical_name(r2[0]):
                raise RuntimeError(
                    "interleaved FASTQ pair names differ: "
                    f"{r1[0].decode(errors='replace').strip()} vs "
                    f"{r2[0].decode(errors='replace').strip()}"
                )
            r1_out.writelines(r1)
            r2_out.writelines(r2)
            count += 1
    return count


def mapped_primary_names(sam_path: Path) -> set[str]:
    names: set[str] = set()
    with sam_path.open("r", encoding="utf-8") as sam:
        for line in sam:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise RuntimeError("invalid Abismal SAM record")
            flag = int(fields[1])
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            names.add(canonical_sam_name(fields[0]))
    return names


def emit_candidates(
    r1_path: Path,
    r2_path: Path,
    keep: set[str],
    r1_out: gzip.GzipFile,
    r2_out: gzip.GzipFile,
) -> int:
    emitted = 0
    with r1_path.open("rb") as r1_in, r2_path.open("rb") as r2_in:
        while True:
            r1 = read_record(r1_in)
            r2 = read_record(r2_in)
            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None:
                raise RuntimeError("temporary paired FASTQ chunks differ in length")
            name1 = canonical_name(r1[0])
            name2 = canonical_name(r2[0])
            if name1 != name2:
                raise RuntimeError("temporary paired FASTQ chunk names differ")
            if name1.decode("utf-8", errors="strict") in keep:
                r1_out.writelines(r1)
                r2_out.writelines(r2)
                emitted += 1
    return emitted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--candidate-r1", required=True)
    parser.add_argument("--candidate-r2", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--chunk-pairs", type=int, default=1_000_000)
    parser.add_argument("--max-edit-distance", type=float, default=0.20)
    args = parser.parse_args()

    if args.threads < 1:
        raise SystemExit("threads must be positive")
    if args.chunk_pairs < 1:
        raise SystemExit("chunk-pairs must be positive")
    if not 0.0 < args.max_edit_distance <= 1.0:
        raise SystemExit("max-edit-distance must be in (0, 1]")

    work_parent = Path(args.workdir)
    work_parent.mkdir(parents=True, exist_ok=True)
    total_pairs = 0
    candidate_pairs = 0
    chunk_no = 0

    with gzip.open(args.candidate_r1, "wb", compresslevel=1) as r1_candidates, gzip.open(
        args.candidate_r2, "wb", compresslevel=1
    ) as r2_candidates:
        while True:
            with tempfile.TemporaryDirectory(prefix="chunk-", dir=work_parent) as temp_dir:
                temp = Path(temp_dir)
                r1_chunk = temp / "R1.fastq"
                r2_chunk = temp / "R2.fastq"
                sam_path = temp / "abismal.sam"
                stats_path = temp / "abismal.stats.yaml"

                count = write_chunk(
                    sys.stdin.buffer, r1_chunk, r2_chunk, args.chunk_pairs
                )
                if count == 0:
                    break
                chunk_no += 1
                total_pairs += count

                command = [
                    "micromamba",
                    "run",
                    "-n",
                    "abismal",
                    "abismal",
                    "map",
                    "-i",
                    args.index,
                    "-t",
                    str(args.threads),
                    "-a",
                    "-m",
                    str(args.max_edit_distance),
                    "-s",
                    str(stats_path),
                    "-o",
                    str(sam_path),
                    str(r1_chunk),
                    str(r2_chunk),
                ]
                subprocess.run(command, check=True)
                if not sam_path.is_file() or sam_path.stat().st_size == 0:
                    raise RuntimeError("Abismal produced no SAM output")

                keep = mapped_primary_names(sam_path)
                emitted = emit_candidates(
                    r1_chunk, r2_chunk, keep, r1_candidates, r2_candidates
                )
                candidate_pairs += emitted
                print(
                    f"Abismal chunk {chunk_no}: input_pairs={count} "
                    f"candidate_pairs={emitted} total_candidates={candidate_pairs}",
                    file=sys.stderr,
                    flush=True,
                )

    if total_pairs == 0:
        raise SystemExit("cleaned interleaved FASTQ contained zero pairs")
    if candidate_pairs == 0:
        raise SystemExit("Abismal retained zero candidate pairs")
    print(
        f"Abismal chunked screen complete: input_pairs={total_pairs} "
        f"candidate_pairs={candidate_pairs} chunks={chunk_no}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
