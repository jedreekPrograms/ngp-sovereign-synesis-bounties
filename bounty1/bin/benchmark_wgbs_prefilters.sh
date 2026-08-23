#!/usr/bin/env bash
set -Eeuo pipefail

# Fail-closed benchmark for a possible faster WGBS PACE prefilter.
#
# This script DOES NOT change the production pipeline. It compares an
# intentionally permissive BSMAP target screen with the current bwa-meth target
# screen on the exact same fastp-cleaned prefix of a paired WGBS sample.
# BSMAP is allowed to retain extra pairs, but it is not eligible for production
# if it misses even one qname retained by the bwa-meth baseline in the tested
# subset.
#
# Usage:
#   benchmark_wgbs_prefilters.sh SAMPLE R1.fastq.gz R2.fastq.gz \
#       PACE_SCREEN.fa THREADS READ_PAIRS OUT_DIR

if [[ $# -ne 7 ]]; then
  echo "usage: $0 SAMPLE R1 R2 PACE_SCREEN_FASTA THREADS READ_PAIRS OUT_DIR" >&2
  exit 2
fi

sample="$1"
r1="$2"
r2="$3"
ref="$4"
threads="$5"
read_pairs="$6"
out_dir="$7"

for path in "$r1" "$r2" "$ref"; do
  [[ -s "$path" ]] || { echo "missing or empty input: $path" >&2; exit 2; }
done
[[ "$threads" =~ ^[0-9]+$ ]] && (( threads >= 2 )) || {
  echo "THREADS must be an integer >=2" >&2
  exit 2
}
[[ "$read_pairs" =~ ^[0-9]+$ ]] && (( read_pairs >= 1000 )) || {
  echo "READ_PAIRS must be an integer >=1000" >&2
  exit 2
}

mkdir -p "$out_dir"
work="$out_dir/work"
mkdir -p "$work"

clean_r1="$work/${sample}.benchmark.clean.R1.fastq.gz"
clean_r2="$work/${sample}.benchmark.clean.R2.fastq.gz"
bwa_names="$out_dir/${sample}.bwameth.names.txt"
bsmap_names="$out_dir/${sample}.bsmap.names.txt"
missed="$out_dir/${sample}.bsmap-missed-vs-bwameth.txt"
extra="$out_dir/${sample}.bsmap-extra-vs-bwameth.txt"
summary="$out_dir/${sample}.summary.tsv"

fastp_threads="$threads"
(( fastp_threads > 4 )) && fastp_threads=4
bsmap_threads=$((threads / 2))
(( bsmap_threads < 1 )) && bsmap_threads=1

now() { date +%s; }

# Create one identical, small, paired input for both mapping strategies. This is
# deliberately a prefix benchmark rather than a random subsample so it can be
# reproduced without scanning the entire 50+ GB source files first.
start=$(now)
fastp \
  --in1 "$r1" \
  --in2 "$r2" \
  --out1 "$clean_r1" \
  --out2 "$clean_r2" \
  --reads_to_process "$read_pairs" \
  --thread "$fastp_threads" \
  --compression 1 \
  --json "$out_dir/${sample}.fastp.json" \
  --html "$out_dir/${sample}.fastp.html" \
  2> "$out_dir/${sample}.fastp.stderr.log"
fastp_seconds=$(( $(now) - start ))

gzip -t "$clean_r1"
gzip -t "$clean_r2"
clean_pairs=$(gzip -cd "$clean_r1" | awk 'END { if (NR % 4 != 0) exit 3; print NR/4 }')
clean_pairs_r2=$(gzip -cd "$clean_r2" | awk 'END { if (NR % 4 != 0) exit 3; print NR/4 }')
[[ "$clean_pairs" == "$clean_pairs_r2" ]] || {
  echo "fastp subset R1/R2 pair counts differ" >&2
  exit 1
}

# Current production baseline: paired, interleaved bwa-meth against the same
# generous PACE screen reference. Keep the whole template whenever either mate
# maps, exactly matching screen_wgbs_pace_candidates.sh semantics.
start=$(now)
python3 /work/bin/stream_interleaved_fastq.py \
  --r1-source "$clean_r1" \
  --r2-source "$clean_r2" \
  2> "$out_dir/${sample}.bwameth.stream.stderr.log" | \
bwameth.py \
  --threads "$threads" \
  --interleaved \
  --reference "$ref" \
  /dev/stdin \
  2> "$out_dir/${sample}.bwameth.stderr.log" | \
samtools view -F 2304 -e '!(flag.unmap && flag.munmap)' - | \
awk '{q=$1; sub(/\/[12]$/, "", q); print q}' | \
LC_ALL=C sort -u > "$bwa_names"
bwameth_seconds=$(( $(now) - start ))

# Experimental screen: map R1 and R2 independently with an intentionally
# permissive BSMAP configuration and union their qnames. Independent mate
# mapping is deliberate: a true target pair is retained even when only one mate
# maps to the small target reference. False positives are acceptable here
# because every retained pair is still realigned to complete hg19 later.
#
# -3   : 3-nucleotide mapping mode (screening only)
# -n 1 : search all four bisulfite strands
# -s 12: shorter seed than default for sensitivity
# -v .12 / -g 3: permissive mismatch/gap settings
# -r 1 : report one best hit; only qname membership matters
bsmap_one_mate() {
  local fq="$1" label="$2" names="$3"
  micromamba run -n bsmap bsmap \
    -a "$fq" \
    -d "$ref" \
    -p "$bsmap_threads" \
    -3 -n 1 -s 12 -v 0.12 -g 3 -r 1 \
    2> "$out_dir/${sample}.bsmap.${label}.stderr.log" | \
  samtools view -F 4 - | \
  awk '{q=$1; sub(/\/[12]$/, "", q); print q}' | \
  LC_ALL=C sort -u > "$names"
}

start=$(now)
bsmap_one_mate "$clean_r1" R1 "$work/bsmap.R1.names.txt" &
p1=$!
bsmap_one_mate "$clean_r2" R2 "$work/bsmap.R2.names.txt" &
p2=$!
wait "$p1"
wait "$p2"
LC_ALL=C sort -u "$work/bsmap.R1.names.txt" "$work/bsmap.R2.names.txt" > "$bsmap_names"
bsmap_seconds=$(( $(now) - start ))

LC_ALL=C comm -23 "$bwa_names" "$bsmap_names" > "$missed"
LC_ALL=C comm -13 "$bwa_names" "$bsmap_names" > "$extra"

bwa_count=$(wc -l < "$bwa_names")
bsmap_count=$(wc -l < "$bsmap_names")
missed_count=$(wc -l < "$missed")
extra_count=$(wc -l < "$extra")

if (( bsmap_seconds > 0 )); then
  speedup=$(awk -v a="$bwameth_seconds" -v b="$bsmap_seconds" 'BEGIN { printf "%.3f", a/b }')
else
  speedup="inf"
fi

status="PASS_SUBSET_RECALL"
if (( missed_count != 0 )); then
  status="FAIL_RECALL"
fi

{
  printf 'metric\tvalue\n'
  printf 'sample\t%s\n' "$sample"
  printf 'requested_pairs\t%s\n' "$read_pairs"
  printf 'clean_pairs\t%s\n' "$clean_pairs"
  printf 'threads_total\t%s\n' "$threads"
  printf 'fastp_seconds\t%s\n' "$fastp_seconds"
  printf 'bwameth_seconds\t%s\n' "$bwameth_seconds"
  printf 'bsmap_seconds\t%s\n' "$bsmap_seconds"
  printf 'bsmap_speedup_vs_bwameth\t%s\n' "$speedup"
  printf 'bwameth_candidate_qnames\t%s\n' "$bwa_count"
  printf 'bsmap_candidate_qnames\t%s\n' "$bsmap_count"
  printf 'bsmap_missed_vs_bwameth\t%s\n' "$missed_count"
  printf 'bsmap_extra_vs_bwameth\t%s\n' "$extra_count"
  printf 'status\t%s\n' "$status"
} > "$summary"

cat "$summary"

# Fail closed. A benchmark miss is a hard failure and must never silently
# promote the experimental screen into production.
if (( missed_count != 0 )); then
  echo "BSMAP screen missed $missed_count bwa-meth candidate qnames; NOT safe to promote." >&2
  exit 1
fi
