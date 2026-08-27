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
# When complete hg19 + official PACE +/-500 bp windows are mounted by either
# the production cohort (/shard) or local runner (/local), this script derives
# the validated +/-2 kb screen itself. The passed SCREEN_REFERENCE_FASTA
# remains a standalone fallback.
#
# Abismal requires two mate filenames and can open/read them sequentially, so
# paired FIFOs can deadlock with an interleaved producer. The cleaned fastp
# stream is therefore processed in bounded temporary chunks: each chunk is
# mapped, candidate pairs are appended to gzip outputs, and the chunk is
# deleted before the next one. Whole cleaned libraries and whole-screen BAMs
# are never materialised.
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

# Derive the benchmark-validated +/-2 kb context from complete hg19 wherever
# the current runner mounts it. This prevents local execution from silently
# falling back to the older +/-1 kb screen that did not achieve zero misses in
# the full-hg19 truth benchmark.
ref_root=""
if [[ -s /shard/ref/hg19.fa && -s /shard/ref/pace.bed ]]; then
  ref_root="/shard/ref"
elif [[ -s /local/ref/hg19.fa && -s /local/ref/pace.bed ]]; then
  ref_root="/local/ref"
fi

if [[ -n "$ref_root" ]]; then
  ctx_bed="/tmp/${sample_id}.pace-screen-2k.merged.bed"
  ctx_fasta="/tmp/${sample_id}.pace-screen-2k.fa"
  python3 - "${ref_root}/pace.bed" "$ctx_bed" <<'PY'
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
    samtools faidx "${ref_root}/hg19.fa" "${chrom}:$((start + 1))-${end}"
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
cleanup() {
  rm -rf "$workdir"
}
trap cleanup EXIT

set +e
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
python3 "${script_dir}/chunked_abismal_screen.py" \
  --index "$screen_index" \
  --threads "$threads" \
  --candidate-r1 "$candidate_r1" \
  --candidate-r2 "$candidate_r2" \
  --workdir "$workdir" \
  --chunk-pairs 4000000 \
  --max-edit-distance 0.20 \
  2> "${sample_id}.pace-screen.abismal.stderr.log"
feed_status=("${PIPESTATUS[@]}")
set -e

{
  echo "stream_interleaved_fastq=${feed_status[0]}"
  echo "fastp=${feed_status[1]}"
  echo "chunked_abismal_screen=${feed_status[2]}"
} > "$status_file"

for code in "${feed_status[@]}"; do
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
