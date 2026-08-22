#!/usr/bin/env bash
set -euo pipefail

# Stream one or two compressed FASTQ sources directly into fastp while
# calculating the archive MD5 on the exact compressed byte stream. Remote
# inputs are never materialized as raw FASTQ files on disk.
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

stream_source() {
  local src="$1"
  case "$src" in
    http://*|https://*|ftp://*)
      curl --fail --location --retry 5 --retry-delay 5 --connect-timeout 30 "$src"
      ;;
    *)
      cat "$src"
      ;;
  esac
}

verify_md5() {
  local expected="$1"
  local actual_file="$2"
  local label="$3"
  [[ -z "$expected" ]] && return 0
  local actual
  actual="$(tr -d '[:space:]' < "$actual_file")"
  if [[ "$actual" != "$expected" ]]; then
    echo "MD5 mismatch for ${label}: expected ${expected}, got ${actual}" >&2
    return 1
  fi
  echo "MD5 OK: ${label} ${actual}"
}

cleanup() {
  rm -f r1.fastp.fifo r1.md5.fifo r2.fastp.fifo r2.md5.fifo
}
trap cleanup EXIT

mkfifo r1.fastp.fifo r1.md5.fifo
(
  stream_source "$r1_source" | tee r1.md5.fifo > r1.fastp.fifo
) &
r1_stream_pid=$!
(
  md5sum < r1.md5.fifo | awk '{print $1}' > r1.actual.md5
) &
r1_md5_pid=$!

if [[ -n "$r2_source" ]]; then
  mkfifo r2.fastp.fifo r2.md5.fifo
  (
    stream_source "$r2_source" | tee r2.md5.fifo > r2.fastp.fifo
  ) &
  r2_stream_pid=$!
  (
    md5sum < r2.md5.fifo | awk '{print $1}' > r2.actual.md5
  ) &
  r2_md5_pid=$!

  fastp \
    -i r1.fastp.fifo \
    -I r2.fastp.fifo \
    -o "${sample_id}.trimmed.R1.fastq.gz" \
    -O "${sample_id}.trimmed.R2.fastq.gz" \
    --json "${sample_id}.fastp.json" \
    --html "${sample_id}.fastp.html" \
    --thread "$threads"

  wait "$r1_stream_pid" "$r1_md5_pid" "$r2_stream_pid" "$r2_md5_pid"
  verify_md5 "$r2_expected" r2.actual.md5 "${sample_id} R2"
else
  fastp \
    -i r1.fastp.fifo \
    -o "${sample_id}.trimmed.R1.fastq.gz" \
    --json "${sample_id}.fastp.json" \
    --html "${sample_id}.fastp.html" \
    --thread "$threads"

  wait "$r1_stream_pid" "$r1_md5_pid"
fi

verify_md5 "$r1_expected" r1.actual.md5 "${sample_id} R1"
