#!/usr/bin/env bash
set -euo pipefail

# Disk-bounded WGBS path.
#
# R1/R2 may be URLs or local compressed FASTQs. stream_interleaved_fastq.py
# validates the exact compressed source bytes while emitting interleaved FASTQ;
# fastp remains a stream; bwa-meth aligns the complete passing-read stream; only
# MAPQ-filtered alignments overlapping the official DunedinPACE probe windows
# are retained.
#
# IMPORTANT: the expensive name-sort is deliberately performed *after* bwa-meth
# has finished. Running bwa-meth and samtools sort concurrently on a 16-GiB
# GitHub runner can exceed RAM for the doubled bisulfite reference. The first
# pipeline therefore writes only the small target BAM, then sorts/deduplicates
# that subset in separate phases. PIPESTATUS and stderr files are retained for
# diagnosis if any stage fails.
#
# Usage:
#   stream_wgbs_bwameth.sh SAMPLE R1 R2 R1_MD5 R2_MD5 \
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

fastp_threads="$threads"
if (( fastp_threads > 2 )); then
  fastp_threads=2
fi
sort_mem="${WGBS_SORT_MEM:-256M}"
target_bam="${sample_id}.mapq${mapq}.pace-targets.unsorted.bam"
name_bam="${sample_id}.pace-targets.name.bam"
fixmate_bam="${sample_id}.pace-targets.fixmate.bam"
position_bam="${sample_id}.pace-targets.position.bam"
final_bam="${sample_id}.mapq${mapq}.pace-targets.deduplicated.bam"
status_file="${sample_id}.pipeline-status.txt"

# Do not let `set -e` hide which member of this long pipeline failed.
set +e
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
  --thread "$fastp_threads" \
  2> "${sample_id}.fastp.stderr.log" | \
bwameth.py \
  --threads "$threads" \
  --interleaved \
  --reference "$reference" \
  /dev/stdin \
  2> "${sample_id}.bwameth.stderr.log" | \
samtools view -bh -q "$mapq" -L "$target_bed" -o "$target_bam" - \
  2> "${sample_id}.samtools-view.stderr.log"
pipe_status=("${PIPESTATUS[@]}")
set -e

{
  echo "stream_interleaved_fastq=${pipe_status[0]}"
  echo "fastp=${pipe_status[1]}"
  echo "bwameth=${pipe_status[2]}"
  echo "samtools_view=${pipe_status[3]}"
} > "$status_file"

pipeline_failed=0
for code in "${pipe_status[@]}"; do
  if (( code != 0 )); then
    pipeline_failed=1
  fi
done
if (( pipeline_failed )); then
  cat "$status_file" >&2
  echo "WGBS streaming pipeline failed; inspect ${sample_id}.*.log" >&2
  exit 1
fi

samtools quickcheck -v "$target_bam"
if [[ "$(samtools view -c "$target_bam")" -eq 0 ]]; then
  echo "No MAPQ${mapq}+ alignments overlapped DunedinPACE target windows" >&2
  exit 1
fi

# Sort only the retained target subset after bwa-meth has released its large
# genome index from RAM. Cap memory per sort thread explicitly.
samtools sort -n -@ "$threads" -m "$sort_mem" -o "$name_bam" "$target_bam"
samtools fixmate -m "$name_bam" "$fixmate_bam"
samtools sort -@ "$threads" -m "$sort_mem" -o "$position_bam" "$fixmate_bam"
samtools markdup -r -@ "$threads" "$position_bam" "$final_bam"
samtools index "$final_bam"
samtools flagstat "$final_bam" > "${sample_id}.pace-targets.flagstat.txt"
samtools quickcheck -v "$final_bam"

rm -f "$target_bam" "$name_bam" "$fixmate_bam" "$position_bam"
printf '%s\n' "$final_bam"
