#!/usr/bin/env bash
set -euo pipefail

# Fast, conservative WGBS candidate screen for DunedinPACE loci.
#
# IMPORTANT: this is only a candidate prefilter. Final candidate alignment,
# MAPQ>=30 filtering, target restriction, duplicate removal and methylation
# calling still use checksum-verified complete hg19.
#
# Validation evidence (GitHub Actions run 32740878122): on 100k authentic
# WT_rep1 pairs, Abismal 3.3.0 against PACE +/-2 kb with `-a -m 0.20` retained
# every one of 603 full-hg19 MAPQ30 PACE truth pairs (0 missed), while taking
# ~12 s and retaining 22,545 candidate pairs. Wider 5/10 kb variants also had
# zero misses; +/-2 kb is therefore the smallest validated context.
#
# When complete hg19 + official PACE +/-500 bp windows are mounted by the
# production cohort, this script derives the validated +/-2 kb screen itself.
# The passed SCREEN_REFERENCE_FASTA remains a standalone fallback.
#
# The cleaned FASTQ stream is split into named pipes because Abismal requires
# separate mate filenames. Abismal SAM is also passed through a named pipe, so
# neither whole cleaned FASTQs nor a whole screen SAM/BAM are materialised.
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

# Production v3 currently prepares a legacy fallback screen before invoking
# this helper. Derive the benchmark-validated context here from complete hg19
# so deployment cannot silently regress to the older +/-1 kb behavior.
if [[ -s /shard/ref/hg19.fa && -s /shard/ref/pace.bed ]]; then
  ctx_bed="/tmp/${sample_id}.pace-screen-2k.merged.bed"
  ctx_fasta="/tmp/${sample_id}.pace-screen-2k.fa"
  python3 - /shard/ref/pace.bed "$ctx_bed" <<'PY'
from collections import defaultdict
from pathlib import Path
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
# pace.bed already carries +/-500 bp; add 1500 bp to obtain +/-2000 total.
extra = 1500
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
  : > "$ctx_fasta"
  while read -r chrom start end; do
    samtools faidx /shard/ref/hg19.fa "${chrom}:$((start + 1))-${end}"
  done < "$ctx_bed" >> "$ctx_fasta"
  test -s "$ctx_fasta"
  screen_reference="$ctx_fasta"
fi

screen_index="${screen_reference}.abismal.idx"
if [[ ! -s "$screen_index" ]]; then
  micromamba run -n abismal abismal idx -t "$threads" \
    "$screen_reference" "$screen_index"
fi
test -s "$screen_index"

candidate_r1="${sample_id}.candidate.R1.fastq.gz"
candidate_r2="${sample_id}.candidate.R2.fastq.gz"
status_file="${sample_id}.pace-screen.pipeline-status.txt"
workdir="$(mktemp -d "/tmp/${sample_id}.abismal-screen.XXXXXX")"
r1_fifo="${workdir}/clean.R1.fastq"
r2_fifo="${workdir}/clean.R2.fastq"
sam_fifo="${workdir}/abismal.sam"
mkfifo "$r1_fifo" "$r2_fifo" "$sam_fifo"
cleanup() {
  rm -rf "$workdir"
}
trap cleanup EXIT

# Consume the SAM FIFO before starting Abismal, then start Abismal before the
# FASTQ splitter. This ordering prevents named-pipe open deadlocks.
set +e
(
  samtools view -h -F SECONDARY,SUPPLEMENTARY \
    -e '!(flag.unmap && flag.munmap)' "$sam_fifo" \
    2> "${sample_id}.pace-screen.samtools-view.stderr.log" | \
  samtools collate -u -O -@ "$threads" - \
    2> "${sample_id}.pace-screen.samtools-collate.stderr.log" | \
  samtools fastq \
    -@ "$fastp_threads" -n \
    -1 "$candidate_r1" -2 "$candidate_r2" \
    -0 /dev/null -s /dev/null - \
    2> "${sample_id}.pace-screen.samtools-fastq.stderr.log"
) &
extract_pid=$!

(
  micromamba run -n abismal abismal map \
    -i "$screen_index" -t "$threads" -a -m 0.20 \
    -s "${sample_id}.pace-screen.abismal.stats.yaml" \
    -o "$sam_fifo" "$r1_fifo" "$r2_fifo" \
    2> "${sample_id}.pace-screen.abismal.stderr.log"
) &
abismal_pid=$!

python3 "${script_dir}/stream_interleaved_fastq.py" \
  --r1-source "$r1_source" --r2-source "$r2_source" \
  --r1-md5 "$r1_md5" --r2-md5 "$r2_md5" \
  2> "${sample_id}.pace-screen.source-stream.log" | \
fastp \
  --stdin --interleaved_in --stdout \
  --json "${sample_id}.pace-screen.fastp.json" \
  --html "${sample_id}.pace-screen.fastp.html" \
  --thread "$fastp_threads" \
  2> "${sample_id}.pace-screen.fastp.stderr.log" | \
python3 "${script_dir}/split_interleaved_fastq.py" \
  --r1-output "$r1_fifo" --r2-output "$r2_fifo" \
  2> "${sample_id}.pace-screen.split.stderr.log"
feed_status=("${PIPESTATUS[@]}")

wait "$abismal_pid"; abismal_status=$?
wait "$extract_pid"; extract_status=$?
set -e

{
  echo "stream_interleaved_fastq=${feed_status[0]}"
  echo "fastp=${feed_status[1]}"
  echo "split_interleaved_fastq=${feed_status[2]}"
  echo "abismal_screen=${abismal_status}"
  echo "samtools_extract_pipeline=${extract_status}"
} > "$status_file"

for code in "${feed_status[@]}" "$abismal_status" "$extract_status"; do
  if (( code != 0 )); then
    cat "$status_file" >&2
    echo "PACE candidate screen failed; inspect ${sample_id}.pace-screen.*.log" >&2
    exit 1
  fi
done

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
printf '%s\n%s\n' "$candidate_r1" "$candidate_r2"
