#!/usr/bin/env python3
"""Stream gzip FASTQ to stdout while validating compressed-source MD5.

Paired reads are emitted interleaved for ``fastp --stdin --interleaved_in``.
Remote HTTP(S) inputs use curl-backed byte-range resume so repeated premature
EOFs from CNCB do not force raw FASTQ staging or restart the compressed stream.
"""

import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from contextlib import ExitStack
from pathlib import Path

REMOTE_SCHEMES = {"http", "https"}


class CurlResumableHTTPReader:
    """Sequential HTTP reader that restarts curl at the next compressed byte."""

    def __init__(self, source, timeout=120, max_retries=2000, retry_delay=1):
        self.source = source
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.offset = 0
        self.total_size = self._probe_total_size()
        self.retries = 0
        self.proc = None
        self._start(0)

    def _probe_total_size(self):
        request = urllib.request.Request(
            self.source,
            headers={
                "User-Agent": "bounty1-fastq-stream/3.0",
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

    def _curl_command(self, offset):
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
            f"{offset}-",
            "--user-agent",
            "bounty1-fastq-stream/3.0",
            self.source,
        ]

    def _finish_process(self):
        if self.proc is None:
            return 0, ""
        if self.proc.stdout is not None:
            self.proc.stdout.close()
        stderr = b""
        if self.proc.stderr is not None:
            stderr = self.proc.stderr.read()
            self.proc.stderr.close()
        rc = self.proc.wait()
        self.proc = None
        return rc, stderr.decode("utf-8", errors="replace").strip()

    def _start(self, offset):
        if self.proc is not None:
            self.proc.kill()
            self._finish_process()
        self.proc = subprocess.Popen(  # noqa: S603
            self._curl_command(offset),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if self.proc.stdout is None:
            raise RuntimeError("curl stdout pipe was not created")

    def _resume(self, reason):
        self.retries += 1
        if self.retries > self.max_retries:
            raise IOError(
                f"exhausted {self.max_retries} curl resume attempts for "
                f"{self.source} at byte {self.offset}"
            ) from reason
        print(
            f"curl resume {self.retries}/{self.max_retries}: {self.source} "
            f"at compressed byte {self.offset}: {reason}",
            file=sys.stderr,
        )
        time.sleep(self.retry_delay)
        self._start(self.offset)

    def read(self, size=-1):
        while True:
            if self.proc is None or self.proc.stdout is None:
                raise RuntimeError("curl reader is closed")
            data = self.proc.stdout.read(size)
            if data:
                self.offset += len(data)
                return data

            rc, stderr = self._finish_process()
            if self.offset >= self.total_size:
                if self.offset != self.total_size:
                    raise IOError(
                        "HTTP stream exceeded expected size: "
                        f"{self.offset} > {self.total_size}"
                    )
                if rc != 0:
                    raise IOError(
                        f"curl exited {rc} after complete object: {stderr or 'no stderr'}"
                    )
                return b""

            self._resume(
                EOFError(
                    f"premature curl EOF at byte {self.offset} of "
                    f"{self.total_size}; rc={rc}; {stderr or 'no stderr'}"
                )
            )

    def close(self):
        if self.proc is not None:
            self.proc.kill()
            self._finish_process()

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
        return CurlResumableHTTPReader(source, timeout=timeout)
    if parsed.scheme == "ftp":
        request = urllib.request.Request(
            source,
            headers={"User-Agent": "bounty1-fastq-stream/3.0"},
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
