#!/usr/bin/env python3
"""Stream gzip FASTQ to stdout while validating compressed-source MD5.

Paired reads are emitted interleaved for ``fastp --stdin --interleaved_in``.
Remote HTTP(S) inputs use bounded parallel curl byte ranges so intermittent
CNCB connection truncation does not force raw FASTQ staging or a restart from
byte zero. Ranges are reassembled strictly in source order and the complete
compressed object is still validated against its published MD5 before success.
"""

import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from pathlib import Path

REMOTE_SCHEMES = {"http", "https"}
RANGE_CHUNK_BYTES = 16 * 1024 * 1024
RANGE_PREFETCH = 3
RANGE_RETRIES = 30


class ParallelCurlRangeReader:
    """Ordered, bounded-memory HTTP reader backed by parallel curl ranges."""

    def __init__(
        self,
        source,
        timeout=120,
        chunk_size=RANGE_CHUNK_BYTES,
        prefetch=RANGE_PREFETCH,
        max_retries=RANGE_RETRIES,
        retry_delay=1,
    ):
        self.source = source
        self.timeout = timeout
        self.chunk_size = chunk_size
        self.prefetch = prefetch
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.total_size = self._probe_total_size()
        self.chunk_count = (self.total_size + self.chunk_size - 1) // self.chunk_size
        self.executor = ThreadPoolExecutor(max_workers=self.prefetch)
        self.pending = {}
        self.next_submit = 0
        self.next_consume = 0
        self.buffer = bytearray()
        self.offset = 0
        self.closed = False
        self._fill_prefetch()

    def _probe_total_size(self):
        request = urllib.request.Request(
            self.source,
            headers={
                "User-Agent": "bounty1-fastq-stream/4.0",
                "Range": "bytes=0-0",
            },
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:  # nosec B310
            status = getattr(response, "status", None)
            content_range = response.headers.get("Content-Range", "")
            if status != 206 or not content_range.startswith("bytes 0-0/"):
                raise IOError(
                    "source does not provide fail-closed HTTP byte ranges: "
                    f"status={status} content-range={content_range!r}"
                )
            total = content_range.rsplit("/", 1)[1]
            if not total.isdigit() or int(total) <= 0:
                raise IOError(f"invalid HTTP object size: {content_range!r}")
            return int(total)

    def _range_bounds(self, index):
        start = index * self.chunk_size
        end = min(self.total_size - 1, start + self.chunk_size - 1)
        return start, end

    def _curl_command(self, start, end):
        connect_timeout = max(5, min(30, self.timeout))
        return [
            "curl",
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "--http1.1",
            "--connect-timeout",
            str(connect_timeout),
            "--speed-time",
            "120",
            "--speed-limit",
            "1024",
            "--range",
            f"{start}-{end}",
            "--user-agent",
            "bounty1-fastq-stream/4.0",
            self.source,
        ]

    def _fetch_chunk(self, index):
        start, end = self._range_bounds(index)
        expected = end - start + 1
        last_error = None
        for attempt in range(1, self.max_retries + 1):
            try:
                completed = subprocess.run(  # noqa: S603
                    self._curl_command(start, end),
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=max(300, self.timeout * 3),
                )
                if completed.returncode != 0:
                    raise IOError(
                        f"curl exit={completed.returncode}: "
                        f"{completed.stderr.decode('utf-8', errors='replace').strip()}"
                    )
                if len(completed.stdout) != expected:
                    raise IOError(
                        f"range length mismatch {start}-{end}: "
                        f"expected {expected}, got {len(completed.stdout)}"
                    )
                return completed.stdout
            except (OSError, subprocess.TimeoutExpired) as exc:
                last_error = exc
                print(
                    f"range retry {attempt}/{self.max_retries}: {self.source} "
                    f"bytes={start}-{end}: {exc}",
                    file=sys.stderr,
                )
                if attempt < self.max_retries:
                    time.sleep(self.retry_delay)
        raise IOError(
            f"failed HTTP range {start}-{end} after {self.max_retries} attempts"
        ) from last_error

    def _fill_prefetch(self):
        while len(self.pending) < self.prefetch and self.next_submit < self.chunk_count:
            index = self.next_submit
            self.pending[index] = self.executor.submit(self._fetch_chunk, index)
            self.next_submit += 1

    def _consume_next_chunk(self):
        if self.next_consume >= self.chunk_count:
            return False
        future = self.pending.pop(self.next_consume, None)
        if future is None:
            raise RuntimeError(f"missing prefetched HTTP range {self.next_consume}")
        data = future.result()
        self.buffer.extend(data)
        self.next_consume += 1
        self._fill_prefetch()
        return True

    def read(self, size=-1):
        if self.closed:
            raise RuntimeError("HTTP range reader is closed")
        if size is None or size < 0:
            while self._consume_next_chunk():
                pass
            data = bytes(self.buffer)
            self.buffer.clear()
            self.offset += len(data)
            return data
        if size == 0:
            return b""

        while len(self.buffer) < size and self.next_consume < self.chunk_count:
            self._consume_next_chunk()

        if not self.buffer:
            if self.offset != self.total_size:
                raise IOError(
                    f"HTTP stream ended at {self.offset} of {self.total_size} compressed bytes"
                )
            return b""

        take = min(size, len(self.buffer))
        data = bytes(self.buffer[:take])
        del self.buffer[:take]
        self.offset += len(data)
        if self.offset > self.total_size:
            raise IOError(
                f"HTTP stream exceeded expected size: {self.offset} > {self.total_size}"
            )
        return data

    def close(self):
        if self.closed:
            return
        self.closed = True
        for future in self.pending.values():
            future.cancel()
        self.pending.clear()
        self.executor.shutdown(wait=False, cancel_futures=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


class HashingReader:
    def __init__(self, raw):
        self.raw = raw
        self.digest = hashlib.md5(usedforsecurity=False)
        self.bytes_read = 0

    def read(self, size=-1):
        data = self.raw.read(size)
        if data:
            self.digest.update(data)
            self.bytes_read += len(data)
        return data

    def drain(self):
        while self.read(1024 * 1024):
            pass

    def hexdigest(self):
        return self.digest.hexdigest()

    def close(self):
        self.raw.close()


def open_source(source, timeout):
    parsed = urllib.parse.urlparse(source)
    if parsed.scheme in REMOTE_SCHEMES:
        if shutil.which("curl") is None:
            raise RuntimeError("curl is required for resumable HTTP FASTQ streaming")
        return ParallelCurlRangeReader(source, timeout=timeout)
    if parsed.scheme == "ftp":
        request = urllib.request.Request(
            source,
            headers={"User-Agent": "bounty1-fastq-stream/4.0"},
        )
        return urllib.request.urlopen(request, timeout=timeout)  # nosec B310
    if parsed.scheme:
        raise ValueError(f"unsupported source scheme: {parsed.scheme}")
    return Path(source).open("rb")


def normalise_read_name(header):
    token = header.strip().split(maxsplit=1)[0]
    if token.startswith(b"@"):
        token = token[1:]
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    return token


def read_record(handle, label):
    header = handle.readline()
    if header == b"":
        return None
    record = [header] + [handle.readline() for _ in range(3)]
    if any(line == b"" for line in record):
        raise EOFError(f"{label}: truncated FASTQ record")
    if not record[0].startswith(b"@"):
        raise ValueError(f"{label}: malformed FASTQ header")
    if not record[2].startswith(b"+"):
        raise ValueError(f"{label}: malformed FASTQ plus line")
    if len(record[1].rstrip()) != len(record[3].rstrip()):
        raise ValueError(f"{label}: sequence and quality lengths differ")
    return record


def verify_md5(label, actual, expected):
    expected = (expected or "").strip().lower()
    if expected and actual.lower() != expected:
        raise ValueError(f"MD5 mismatch for {label}: expected {expected}, got {actual}")
    if expected:
        print(f"MD5 OK: {label} {actual}", file=sys.stderr)


def stream_fastq(r1_source, r2_source, out, r1_md5="", r2_md5="", timeout=120):
    with ExitStack() as stack:
        raw1 = stack.enter_context(open_source(r1_source, timeout))
        hashed1 = HashingReader(raw1)
        gz1 = gzip.GzipFile(fileobj=hashed1, mode="rb")
        stack.callback(gz1.close)

        hashed2 = None
        gz2 = None
        if r2_source:
            raw2 = stack.enter_context(open_source(r2_source, timeout))
            hashed2 = HashingReader(raw2)
            gz2 = gzip.GzipFile(fileobj=hashed2, mode="rb")
            stack.callback(gz2.close)

        record_index = 0
        while True:
            r1 = read_record(gz1, f"R1 record {record_index + 1}")
            r2 = read_record(gz2, f"R2 record {record_index + 1}") if gz2 else None
            if r1 is None:
                if r2 is not None:
                    raise ValueError("R2 contains more records than R1")
                break
            if gz2 and r2 is None:
                raise ValueError("R1 contains more records than R2")
            if r2 is not None and normalise_read_name(r1[0]) != normalise_read_name(r2[0]):
                raise ValueError(
                    f"pair-name mismatch at record {record_index + 1}: "
                    f"{r1[0].strip()!r} != {r2[0].strip()!r}"
                )
            out.writelines(r1)
            if r2 is not None:
                out.writelines(r2)
            record_index += 1

        hashed1.drain()
        if hashed2 is not None:
            hashed2.drain()

        verify_md5("R1", hashed1.hexdigest(), r1_md5)
        if hashed2 is not None:
            verify_md5("R2", hashed2.hexdigest(), r2_md5)

    return record_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1-source", required=True)
    parser.add_argument("--r2-source", default="")
    parser.add_argument("--r1-md5", default="")
    parser.add_argument("--r2-md5", default="")
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()

    count = stream_fastq(
        args.r1_source,
        args.r2_source,
        sys.stdout.buffer,
        r1_md5=args.r1_md5,
        r2_md5=args.r2_md5,
        timeout=args.timeout,
    )
    sys.stdout.buffer.flush()
    print(f"Streamed FASTQ records: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
