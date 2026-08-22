#!/usr/bin/env python3
"""Extract the first N paired FASTQ records from public gzip URLs.

This is for CI/technical smoke tests only. It deliberately reads only a prefix
of each source FASTQ, so archive MD5 values for the complete files are not
expected to match. Production analysis uses the full streaming path and checks
source MD5 values.
"""
import argparse
import gzip
import hashlib
import json
import urllib.parse
import urllib.request
from contextlib import ExitStack
from pathlib import Path


def normalise_read_name(header: bytes) -> bytes:
    token = header.strip().split(maxsplit=1)[0]
    if token.startswith(b'@'):
        token = token[1:]
    if token.endswith(b'/1') or token.endswith(b'/2'):
        token = token[:-2]
    return token


def read_record(handle, label):
    record = [handle.readline() for _ in range(4)]
    if any(line == b'' for line in record):
        raise EOFError(f'{label}: source ended before a complete FASTQ record')
    if not record[0].startswith(b'@'):
        raise ValueError(f'{label}: malformed FASTQ header: {record[0][:80]!r}')
    if not record[2].startswith(b'+'):
        raise ValueError(f'{label}: malformed FASTQ plus line')
    if len(record[1].rstrip()) != len(record[3].rstrip()):
        raise ValueError(f'{label}: sequence and quality lengths differ')
    return record


def open_gzip_source(url, timeout):
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme == 'https':
        request = urllib.request.Request(
            url,
            headers={'User-Agent': 'bounty1-real-smoke/1.0'},
        )
        response = urllib.request.urlopen(request, timeout=timeout)  # nosec B310
    elif parsed.scheme == 'file':
        local_path = Path(urllib.request.url2pathname(parsed.path))
        response = local_path.open('rb')
    else:
        raise ValueError(f'unsupported FASTQ URL scheme: {parsed.scheme!r}')
    return response, gzip.GzipFile(fileobj=response, mode='rb')


def extract(url1, url2, out1, out2, reads, timeout):
    if reads < 1:
        raise ValueError('--reads must be >= 1')

    with ExitStack() as stack:
        response1, in1 = open_gzip_source(url1, timeout)
        stack.callback(response1.close)
        stack.callback(in1.close)
        response2, in2 = open_gzip_source(url2, timeout)
        stack.callback(response2.close)
        stack.callback(in2.close)

        raw1 = stack.enter_context(Path(out1).open('wb'))
        raw2 = stack.enter_context(Path(out2).open('wb'))
        gz1 = gzip.GzipFile(fileobj=raw1, mode='wb', mtime=0)
        gz2 = gzip.GzipFile(fileobj=raw2, mode='wb', mtime=0)
        stack.callback(gz1.close)
        stack.callback(gz2.close)

        for index in range(reads):
            r1 = read_record(in1, f'R1 record {index + 1}')
            r2 = read_record(in2, f'R2 record {index + 1}')
            if normalise_read_name(r1[0]) != normalise_read_name(r2[0]):
                raise ValueError(
                    f'pair-name mismatch at record {index + 1}: '
                    f'{r1[0].strip()!r} != {r2[0].strip()!r}'
                )
            for line in r1:
                gz1.write(line)
            for line in r2:
                gz2.write(line)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--url1', required=True)
    parser.add_argument('--url2', required=True)
    parser.add_argument('--out1', required=True)
    parser.add_argument('--out2', required=True)
    parser.add_argument('--reads', type=int, default=250_000)
    parser.add_argument('--timeout', type=int, default=120)
    parser.add_argument('--metadata', default='')
    args = parser.parse_args()

    extract(args.url1, args.url2, args.out1, args.out2, args.reads, args.timeout)
    metadata = {
        'source_url_1': args.url1,
        'source_url_2': args.url2,
        'records_per_mate': args.reads,
        'partial_source': True,
        'full_source_md5_checked': False,
        'output_sha256_1': sha256(args.out1),
        'output_sha256_2': sha256(args.out2),
    }
    if args.metadata:
        Path(args.metadata).write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
