#!/usr/bin/env bash
set -euo pipefail

# Fast, conservative WGBS candidate screen for DunedinPACE loci.
#
# This is NOT the final alignment. The complete paired FASTQ source is first
# validated against its published MD5 and passed once through fastp. Cleaned
# reads are aligned to a deliberately generous target-only bisulfite reference.
# Any primary pair for which either mate maps to that screen reference is kept.
# The resulting candidate FASTQs must subsequently be aligned against the full
# hg19 reference; final MAPQ filtering, target restriction, duplicate removal,
# and methylation calling therefore still use the complete genome reference.
#
# In the production cohort /shard/ref/pace.bed contains the official PACE
# windows at +/-500 bp. We expand those windows by another 500 bp before the
# screen, giving +/-1000 bp total around each probe. This is conservative for
# 150-bp WGBS reads while making the screen reference far smaller than the
# historical +/-10 kb fallback reference. The supplied SCREEN_REFERENCE_FASTA
# remains the fallback for standalone use.
#
# Usage:
#   screen_wgbs_pace_candidates.sh SAMPLE R1 R2 R1_MD5 R2_MD5 \
#       SCREEN_REFERENCE_FASTA THREADS

if [[ $# -ne 7 ]]; then
  echo "usage: $0 SAMPLE R1 R2 R1_MD5 R2_MD5 SCREEN_REF_FASTA THREADS" >&2
  exit 2
fi

sample_id="$1"
r1_source="$2"
r2_source="$3"
r1_md5="${4,,}"
r2_md5="${5,,}"
screen_reference="$6"
threads="$7"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "$r2_source" ]]; then
  echo "PACE candidate screen requires paired-end WGBS" >&2
  exit 2
fi
if [[ ! -s "$screen_reference" ]]; then
  echo "screen reference FASTA not found or empty: $screen_reference" >&2
  exit 2
fi
if ! [[ "$threads" =~ ^[0-9]+$ ]] || (( threads < 1 )); then
  echo "THREADS must be a positive integer" >&2
  exit 2
fi

fastp_threads="$threads"
if (( fastp_threads > 2 )); then
  fastp_threads=2
fi

# Build a much smaller but still conservative production screening reference
# whenever the full hg19 reference and the final +/-500 PACE windows are
# available. A read overlapping the final window can extend by at most its read
# length beyond that window; the additional 500 bp margin is therefore safely
# larger than the 150-bp source reads used by HRA003336.
if [[ -s /shard/ref/hg19.fa && -s /shard/ref/pace.bed ]]; then
  narrow_bed="/tmp/${sample_id}.pace-screen-1000.merged.bed"
  narrow_reference="/tmp/${sample_id}.pace-screen-1000.fa"
  python3 - /shard/ref/pace.bed "$narrow_bed" <<'PY'
from collections import defaultdict
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
extra = 500
by_chrom = defaultdict(list)
for line in src.read_text(encoding="utf-8").splitlines():
    if not line or line.startswith("#"):
        continue
    fields = line.split("\t")
    chrom = fields[0]
    start = max(0, int(fields[1]) - extra)
    end = int(fields[2]) + extra
    by_chrom[chrom].append((start, end))

with dst.open("w", encoding="utf-8") as out:
    for chrom in sorted(by_chrom):
        merged = []
        for start, end in sorted(by_chrom[chrom]):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        for start, end in merged:
            out.write(f"{chrom}\t{start}\t{end}\n")
PY

  : > "$narrow_reference"
  while read -r chrom start end; do
    samtools faidx /shard/ref/hg19.fa "${chrom}:$((start + 1))-${end}"
  done < "$narrow_bed" >> "$narrow_reference"
  test -s "$narrow_reference"
  samtools faidx "$narrow_reference"
  bwameth.py index "$narrow_reference"
  screen_reference="$narrow_reference"
fi

# Keep the large transient screen BAM outside /shard so a timeout/cancellation
# cannot accidentally upload gigabytes of a partial BAM as diagnostic evidence.
screen_bam="/tmp/${sample_id}.pace-screen.primary-any-mapped.bam"
candidate_r1="${sample_id}.candidate.R1.fastq.gz"
candidate_r2="${sample_id}.candidate.R2.fastq.gz"
status_file="${sample_id}.pace-screen.pipeline-status.txt"
trap 'rm -f "$screen_bam"' EXIT INT TERM

# `flag.unmap && flag.munmap` is true only when both ends of the template are
# unmapped. Negating it keeps the complete primary pair whenever either mate
# maps to the generous screen reference. Secondary/supplementary records are
# excluded before reconstructing FASTQ.
set +e
python3 "${script_dir}/stream_interleaved_fastq.py" \
  --r1-source "$r1_source" \
  --r2-source "$r2_source" \
  --r1-md5 "$r1_md5" \
  --r2-md5 "$r2_md5" \
  2> "${sample_id}.pace-screen.source-stream.log" | \
fastp \
  --stdin \
  --interleaved_in \
  --stdout \
  --json "${sample_id}.pace-screen.fastp.json" \
  --html "${sample_id}.pace-screen.fastp.html" \
  --thread "$fastp_threads" \
  2> "${sample_id}.pace-screen.fastp.stderr.log" | \
bwameth.py \
  --threads "$threads" \
  --interleaved \
  --reference "$screen_reference" \
  /dev/stdin \
  2> "${sample_id}.pace-screen.bwameth.stderr.log" | \
samtools view -bh \
  -F SECONDARY,SUPPLEMENTARY \
  -e '!(flag.unmap && flag.munmap)' \
  -o "$screen_bam" - \
  2> "${sample_id}.pace-screen.samtools-view.stderr.log"
pipe_status=("${PIPESTATUS[@]}")
set -e

{
  echo "stream_interleaved_fastq=${pipe_status[0]}"
  echo "fastp=${pipe_status[1]}"
  echo "bwameth_screen=${pipe_status[2]}"
  echo "samtools_view=${pipe_status[3]}"
} > "$status_file"

for code in "${pipe_status[@]}"; do
  if (( code != 0 )); then
    cat "$status_file" >&2
    echo "PACE candidate screen failed; inspect ${sample_id}.pace-screen.*.log" >&2
    exit 1
  fi
done

samtools quickcheck -v "$screen_bam"
if [[ "$(samtools view -c "$screen_bam")" -eq 0 ]]; then
  echo "PACE candidate screen retained zero primary alignments" >&2
  exit 1
fi

# Name-collate before FASTQ reconstruction. `samtools fastq` restores reads to
# FASTQ orientation, so the later full-hg19 alignment receives ordinary paired
# reads rather than target-reference-oriented sequences.
samtools collate -u -O -@ "$threads" "$screen_bam" | \
samtools fastq \
  -@ "$fastp_threads" \
  -n \
  -1 "$candidate_r1" \
  -2 "$candidate_r2" \
  -0 /dev/null \
  -s /dev/null \
  - 2> "${sample_id}.pace-screen.samtools-fastq.stderr.log"

gzip -t "$candidate_r1"
gzip -t "$candidate_r2"

r1_lines=$(gzip -cd "$candidate_r1" | wc -l)
r2_lines=$(gzip -cd "$candidate_r2" | wc -l)
if (( r1_lines == 0 || r2_lines == 0 || r1_lines % 4 != 0 || r2_lines % 4 != 0 )); then
  echo "candidate FASTQ record structure is invalid" >&2
  exit 1
fi
if (( r1_lines != r2_lines )); then
  echo "candidate R1/R2 record counts differ" >&2
  exit 1
fi

echo $((r1_lines / 4)) > "${sample_id}.candidate-pair-count.txt"
du -h "$candidate_r1" "$candidate_r2" > "${sample_id}.candidate-fastq-size.txt"

rm -f "$screen_bam"
printf '%s\n%s\n' "$candidate_r1" "$candidate_r2"
