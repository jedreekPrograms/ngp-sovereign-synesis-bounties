#!/usr/bin/env python3
"""Stream gzip FASTQ to stdout while validating source MD5.

For paired-end input, records are emitted as an interleaved FASTQ stream
(R1 record, R2 record, ...), which can be consumed by
`fastp --stdin --interleaved_in` without named pipes or raw-file staging.
MD5 is calculated over the exact compressed source bytes, not decompressed
FASTQ content.
"""
import argparse
import gzip
import hashlib
import sys
import urllib.parse
import urllib.request
from contextlib import ExitStack
from pathlib import Path

REMOTE_SCHEMES = {'http', 'https', 'ftp'}


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

    def tell(self):
        return self.bytes_read

    def drain(self):
        while True:
            data = self.read(1024 * 1024)
            if not data:
                break

    def hexdigest(self):
        return self.digest.hexdigest()

    def close(self):
        self.raw.close()


def open_source(source, timeout):
    parsed = urllib.parse.urlparse(source)
    if parsed.scheme in REMOTE_SCHEMES:
        request = urllib.request.Request(
            source,
            headers={'User-Agent': 'bounty1-fastq-stream/1.0'},
        )
        return urllib.request.urlopen(request, timeout=timeout)  # nosec B310
    if parsed.scheme:
        raise ValueError(f'unsupported source scheme: {parsed.scheme}')
    return Path(source).open('rb')


def normalise_read_name(header):
    token = header.strip().split(maxsplit=1)[0]
    if token.startswith(b'@'):
        token = token[1:]
    if token.endswith((b'/1', b'/2')):
        token = token[:-2]
    return token


def read_record(handle, label):
    header = handle.readline()
    if header == b'':
        return None
    record = [header] + [handle.readline() for _ in range(3)]
    if any(line == b'' for line in record):
        raise EOFError(f'{label}: truncated FASTQ record')
    if not record[0].startswith(b'@'):
        raise ValueError(f'{label}: malformed FASTQ header')
    if not record[2].startswith(b'+'):
        raise ValueError(f'{label}: malformed FASTQ plus line')
    if len(record[1].rstrip()) != len(record[3].rstrip()):
        raise ValueError(f'{label}: sequence and quality lengths differ')
    return record


def verify_md5(label, actual, expected):
    expected = (expected or '').strip().lower()
    if expected and actual.lower() != expected:
        raise ValueError(f'MD5 mismatch for {label}: expected {expected}, got {actual}')
    if expected:
        print(f'MD5 OK: {label} {actual}', file=sys.stderr)


def stream_fastq(r1_source, r2_source, out, r1_md5='', r2_md5='', timeout=120):
    with ExitStack() as stack:
        raw1 = stack.enter_context(open_source(r1_source, timeout))
        hashed1 = HashingReader(raw1)
        gz1 = gzip.GzipFile(fileobj=hashed1, mode='rb')
        stack.callback(gz1.close)

        hashed2 = None
        gz2 = None
        if r2_source:
            raw2 = stack.enter_context(open_source(r2_source, timeout))
            hashed2 = HashingReader(raw2)
            gz2 = gzip.GzipFile(fileobj=hashed2, mode='rb')
            stack.callback(gz2.close)

        record_index = 0
        while True:
            r1 = read_record(gz1, f'R1 record {record_index + 1}')
            r2 = read_record(gz2, f'R2 record {record_index + 1}') if gz2 else None
            if r1 is None:
                if r2 is not None:
                    raise ValueError('R2 contains more records than R1')
                break
            if gz2 and r2 is None:
                raise ValueError('R1 contains more records than R2')
            if r2 is not None and normalise_read_name(r1[0]) != normalise_read_name(r2[0]):
                raise ValueError(
                    f'pair-name mismatch at record {record_index + 1}: '
                    f'{r1[0].strip()!r} != {r2[0].strip()!r}'
                )
            out.writelines(r1)
            if r2 is not None:
                out.writelines(r2)
            record_index += 1

        # Hash any compressed trailer bytes that gzip may not have consumed.
        hashed1.drain()
        if hashed2 is not None:
            hashed2.drain()

        verify_md5('R1', hashed1.hexdigest(), r1_md5)
        if hashed2 is not None:
            verify_md5('R2', hashed2.hexdigest(), r2_md5)

    return record_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r1-source', required=True)
    parser.add_argument('--r2-source', default='')
    parser.add_argument('--r1-md5', default='')
    parser.add_argument('--r2-md5', default='')
    parser.add_argument('--timeout', type=int, default=120)
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
    print(f'Streamed FASTQ records: {count}', file=sys.stderr)


if __name__ == '__main__':
    main()
