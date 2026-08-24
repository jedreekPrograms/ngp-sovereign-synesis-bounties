#!/usr/bin/env python3
"""Split an interleaved FASTQ stream into paired FIFO/file outputs.

This helper deliberately performs no record transformation. It validates basic
FASTQ framing and matching read names, then writes R1/R2 records separately so
mappers that require two filenames can still consume a fully streaming fastp
output without materialising cleaned whole-library FASTQs on disk.
"""

from __future__ import annotations

import argparse
import os
import sys
import threading


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


def open_paired_outputs(r1_path: str, r2_path: str):
    """Open both FIFO writers concurrently to avoid paired-reader deadlocks.

    Some mappers open R1 for reading and then block while opening R2. If this
    process opens R1 for writing first and starts filling it before opening R2,
    both sides can deadlock once the R1 pipe buffer fills. Opening the two write
    ends concurrently lets both FIFO handshakes complete before any reads are
    streamed. Regular files also work with the same code path.
    """

    paths = [r1_path, r2_path]
    fds: list[int | None] = [None, None]
    errors: list[BaseException] = []

    def opener(index: int) -> None:
        try:
            fds[index] = os.open(paths[index], os.O_WRONLY | os.O_CREAT, 0o666)
        except BaseException as exc:  # propagated in the main thread below
            errors.append(exc)

    threads = [threading.Thread(target=opener, args=(i,), daemon=True) for i in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    if errors:
        for fd in fds:
            if fd is not None:
                os.close(fd)
        raise errors[0]
    if any(fd is None for fd in fds):
        raise RuntimeError("failed to open paired FASTQ outputs")

    return os.fdopen(fds[0], "wb", buffering=0), os.fdopen(fds[1], "wb", buffering=0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1-output", required=True)
    parser.add_argument("--r2-output", required=True)
    args = parser.parse_args()

    count = 0
    r1_out, r2_out = open_paired_outputs(args.r1_output, args.r2_output)
    with r1_out, r2_out:
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
