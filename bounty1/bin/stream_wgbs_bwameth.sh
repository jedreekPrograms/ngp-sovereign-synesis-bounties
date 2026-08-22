#!/usr/bin/env bash
set -euo pipefail

# Full-source WGBS streaming path.
#
# The compressed FASTQs are never staged. stream_interleaved_fastq.py validates
# archive MD5 while emitting interleaved FASTQ; fastp QC/trimming remains a
# stream; bwa-meth aligns the complete passing-read stream; only MAPQ-filtered
# alignments overlapping the official DunedinPACE probe windows are retained.
# The retained subset is then duplicate-marked and coordinate sorted for
# MethylDackel. Alignment is therefore genome-wide while disk use is bounded by
# the small target subset needed to calculate DunedinPACE.
#
# Usage:
#   stream_wgbs_bwameth.sh SAMPLE R1_URL R2_URL R1_MD5 R2_MD5 \
#       REFERENCE_FASTA TARGET_WINDOWS_BED MAPQ THREADS

if [[ $# -ne 9 ]]; then
  echo "usage: $0 SAMPLE R1 R2 R1_MD5 R2_MD5 REF_FASTA TARGET_BED MAPQ THREADS" >&2
  exit 2
fi

sample_id="$1"
r1_source="$2"
r2_source="$3"
r1_md5="${4,,}"
r2_md5="${5,,}"
reference="$6"
target_bed="$7"
mapq="$8"
threads="$9"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "$r2_source" ]]; then
  echo "streaming bwa-meth path currently requires paired-end WGBS" >&2
  exit 2
fi
if [[ ! -s "$reference" ]]; then
  echo "reference FASTA not found or empty: $reference" >&2
  exit 2
fi
if [[ ! -s "$target_bed" ]]; then
  echo "target BED not found or empty: $target_bed" >&2
  exit 2
fi
if ! [[ "$mapq" =~ ^[0-9]+$ ]] || (( mapq < 0 )); then
  echo "MAPQ must be a non-negative integer" >&2
  exit 2
fi
if ! [[ "$threads" =~ ^[0-9]+$ ]] || (( threads < 1 )); then
  echo "THREADS must be a positive integer" >&2
  exit 2
fi

name_bam="${sample_id}.pace-targets.name.bam"
fixmate_bam="${sample_id}.pace-targets.fixmate.bam"
position_bam="${sample_id}.pace-targets.position.bam"
final_bam="${sample_id}.mapq${mapq}.pace-targets.deduplicated.bam"

python3 "${script_dir}/stream_interleaved_fastq.py" \
  --r1-source "$r1_source" \
  --r2-source "$r2_source" \
  --r1-md5 "$r1_md5" \
  --r2-md5 "$r2_md5" \
  2> "${sample_id}.source-stream.log" | \
fastp \
  --stdin \
  --interleaved_in \
  --stdout \
  --json "${sample_id}.fastp.json" \
  --html "${sample_id}.fastp.html" \
  --thread "$threads" \
  2> "${sample_id}.fastp.stderr.log" | \
bwameth.py \
  --threads "$threads" \
  --interleaved \
  --reference "$reference" \
  /dev/stdin \
  2> "${sample_id}.bwameth.stderr.log" | \
samtools view -bh -q "$mapq" -L "$target_bed" - | \
samtools sort -n -@ "$threads" -o "$name_bam"

samtools quickcheck -v "$name_bam"
if [[ "$(samtools view -c "$name_bam")" -eq 0 ]]; then
  echo "No MAPQ${mapq}+ alignments overlapped DunedinPACE target windows" >&2
  exit 1
fi

samtools fixmate -m "$name_bam" "$fixmate_bam"
samtools sort -@ "$threads" -o "$position_bam" "$fixmate_bam"
samtools markdup -r -@ "$threads" "$position_bam" "$final_bam"
samtools index "$final_bam"
samtools flagstat "$final_bam" > "${sample_id}.pace-targets.flagstat.txt"
samtools quickcheck -v "$final_bam"

rm -f "$name_bam" "$fixmate_bam" "$position_bam"
printf '%s\n' "$final_bam"
