#!/usr/bin/env bash
set -euo pipefail

# Stream one or two gzip-compressed FASTQ sources through fastp while
# validating archive MD5 on the exact compressed bytes. Paired input is
# converted to an interleaved plain-FASTQ stream and read by fastp from STDIN,
# avoiding named-pipe seek/deadlock behaviour and avoiding raw FASTQ staging.
#
# Usage:
#   stream_fastp.sh SAMPLE_ID R1_SOURCE R2_SOURCE R1_MD5 R2_MD5 THREADS

if [[ $# -ne 6 ]]; then
  echo "usage: $0 SAMPLE_ID R1_SOURCE R2_SOURCE R1_MD5 R2_MD5 THREADS" >&2
  exit 2
fi

sample_id="$1"
r1_source="$2"
r2_source="$3"
r1_expected="${4,,}"
r2_expected="${5,,}"
threads="$6"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "$r2_source" ]]; then
  python3 "${script_dir}/stream_interleaved_fastq.py" \
    --r1-source "$r1_source" \
    --r2-source "$r2_source" \
    --r1-md5 "$r1_expected" \
    --r2-md5 "$r2_expected" | \
  fastp \
    --stdin \
    --interleaved_in \
    -o "${sample_id}.trimmed.R1.fastq.gz" \
    -O "${sample_id}.trimmed.R2.fastq.gz" \
    --json "${sample_id}.fastp.json" \
    --html "${sample_id}.fastp.html" \
    --thread "$threads"
else
  python3 "${script_dir}/stream_interleaved_fastq.py" \
    --r1-source "$r1_source" \
    --r1-md5 "$r1_expected" | \
  fastp \
    --stdin \
    -o "${sample_id}.trimmed.R1.fastq.gz" \
    --json "${sample_id}.fastp.json" \
    --html "${sample_id}.fastp.html" \
    --thread "$threads"
fi
