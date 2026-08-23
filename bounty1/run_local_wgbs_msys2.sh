#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

# Resumable local WGBS runner for Windows + MSYS2 UCRT64 + Docker Desktop.
# Heavy raw files are processed one sample at a time and deleted after a
# validated Bismark-compatible coverage file is produced. Completed samples
# are skipped on rerun and partial FASTQ downloads are resumed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SAMPLESHEET="${SCRIPT_DIR}/resources/samplesheet.csv"

THREADS="${THREADS:-8}"
MIN_FREE_GB="${MIN_FREE_GB:-120}"
SAMPLE_REGEX="${SAMPLE_REGEX:-^(WT|SIRT[1-7])_rep[12]_WGBS$}"
LOCAL_ROOT="${BOUNTY1_LOCAL_ROOT:-${HOME}/bounty1-wgbs}"

REF_DIR="${LOCAL_ROOT}/ref"
WORK_DIR="${LOCAL_ROOT}/work"
RESULTS_DIR="${LOCAL_ROOT}/results"
COV_DIR="${RESULTS_DIR}/cov"
QC_DIR="${RESULTS_DIR}/qc"
LOG_DIR="${LOCAL_ROOT}/logs"
RUN_LOG="${LOG_DIR}/run.log"

mkdir -p "${REF_DIR}" "${WORK_DIR}" "${COV_DIR}" "${QC_DIR}" "${LOG_DIR}"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { printf '[%s] %s\n' "$(timestamp)" "$*" | tee -a "${RUN_LOG}"; }
die() { log "ERROR: $*"; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Brakuje polecenia '$1'."; }

for cmd in git docker curl md5sum gzip awk sort cygpath df find tee seq; do
  require_cmd "$cmd"
done

[[ "${THREADS}" =~ ^[0-9]+$ ]] && (( THREADS >= 1 )) || die "THREADS musi byc dodatnia liczba calkowita."
[[ "${MIN_FREE_GB}" =~ ^[0-9]+$ ]] && (( MIN_FREE_GB >= 50 )) || die "MIN_FREE_GB musi byc >= 50."
[[ -s "${SAMPLESHEET}" ]] || die "Nie znaleziono ${SAMPLESHEET}."

docker info >/dev/null 2>&1 || die "Docker Desktop nie odpowiada. Uruchom Docker Desktop i poczekaj az bedzie gotowy."

# Do not let MSYS2 rewrite Linux paths used inside containers.
export MSYS_NO_PATHCONV=1
export MSYS2_ARG_CONV_EXCL='*'

to_win_path() { cygpath -am "$1"; }
LOCAL_ROOT_WIN="$(to_win_path "${LOCAL_ROOT}")"
BOUNTY_CONTEXT_WIN="$(to_win_path "${SCRIPT_DIR}")"

GIT_SHA="$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || true)"
[[ -n "${GIT_SHA}" ]] || die "Nie moge odczytac aktualnego commita Git."
IMAGE="${BOUNTY1_IMAGE:-bounty1-pace:local-${GIT_SHA:0:12}}"

docker_local() {
  docker run --rm \
    --mount "type=bind,source=${LOCAL_ROOT_WIN},target=/local" \
    "$@"
}

check_disk() {
  local free_kb free_gb
  free_kb="$(df -Pk "${LOCAL_ROOT}" | awk 'NR==2 {print $4}')"
  [[ "${free_kb}" =~ ^[0-9]+$ ]] || die "Nie moge odczytac wolnego miejsca."
  free_gb=$((free_kb / 1024 / 1024))
  log "Wolne miejsce: ${free_gb} GiB (minimum bezpieczenstwa: ${MIN_FREE_GB} GiB)."
  (( free_gb >= MIN_FREE_GB )) || die "Za malo wolnego miejsca; zatrzymuje sie bez kasowania Twoich danych/cache."
}

valid_gzip_md5() {
  local path="$1" expected="${2,,}" actual
  [[ -s "${path}" ]] || return 1
  gzip -t "${path}" >/dev/null 2>&1 || return 1
  actual="$(md5sum "${path}" | awk '{print tolower($1)}')"
  [[ "${actual}" == "${expected}" ]]
}

download_with_resume() {
  local label="$1" url="$2" path="$3" expected="${4,,}"
  local marker="${path}.md5ok" attempt status actual

  if [[ -s "${marker}" ]] && grep -qx "${expected}" "${marker}" && gzip -t "${path}" >/dev/null 2>&1; then
    log "${label}: kompletne dane sa juz na dysku; pomijam download."
    return 0
  fi
  if [[ -s "${path}" ]] && valid_gzip_md5 "${path}" "${expected}"; then
    printf '%s\n' "${expected}" > "${marker}"
    log "${label}: istniejacy plik przeszedl MD5 i gzip; pomijam download."
    return 0
  fi

  for attempt in $(seq 1 100); do
    log "${label}: download/resume proba ${attempt}/100."
    set +e
    curl --fail --location --show-error --http1.1 \
      --retry 5 --retry-delay 5 --retry-all-errors \
      --connect-timeout 30 --speed-time 120 --speed-limit 1024 \
      --continue-at - --output "${path}" "${url}"
    status=$?
    set -e

    if (( status == 0 )); then
      actual="$(md5sum "${path}" | awk '{print tolower($1)}')"
      if [[ "${actual}" == "${expected}" ]] && gzip -t "${path}" >/dev/null 2>&1; then
        printf '%s\n' "${expected}" > "${marker}"
        log "${label}: gotowe; MD5=${actual}."
        return 0
      fi
      log "${label}: transfer zakonczony, ale MD5/gzip jeszcze sie nie zgadza; wznawiam."
    else
      log "${label}: curl exit=${status}; zachowuje czesciowy plik i wznawiam."
    fi
    sleep 5
  done
  return 1
}

build_image() {
  if docker image inspect "${IMAGE}" >/dev/null 2>&1; then
    log "Obraz ${IMAGE} juz istnieje; pomijam build."
    return 0
  fi
  log "Buduje przypiety obraz ${IMAGE}. Pierwszy build moze potrwac."
  docker build --pull -t "${IMAGE}" "${BOUNTY_CONTEXT_WIN}" 2>&1 | tee -a "${RUN_LOG}"
}

prepare_reference() {
  local ready="${REF_DIR}/.hg19-pace-ready-v1" ref_log="${LOG_DIR}/reference.log"
  if [[ -s "${ready}" && -s "${REF_DIR}/hg19.fa" && -s "${REF_DIR}/hg19.fa.fai" \
        && -s "${REF_DIR}/pace.bed" && -s "${REF_DIR}/pace-screen.fa" ]]; then
    log "Referencja hg19 + PACE jest juz gotowa."
    return 0
  fi

  log "Przygotowuje hg19 i indeksy tylko raz. Szczegoly: ${ref_log}"
  rm -f "${ready}"
  {
    set -Eeuo pipefail
    if [[ ! -s "${REF_DIR}/hg19.fa" ]]; then
      curl --fail --location --show-error --retry 20 --retry-delay 5 --retry-all-errors \
        --continue-at - --output "${REF_DIR}/hg19.fa.gz" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
      curl --fail --location --show-error --retry 20 --retry-delay 5 --retry-all-errors \
        --output "${REF_DIR}/md5sum.txt" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt
      (
        cd "${REF_DIR}"
        grep -E '[[:space:]]hg19.fa.gz$' md5sum.txt > hg19.fa.gz.md5
        md5sum --check hg19.fa.gz.md5
        gzip -t hg19.fa.gz
        gzip -dc hg19.fa.gz > hg19.fa.tmp
        mv hg19.fa.tmp hg19.fa
        rm -f hg19.fa.gz
      )
    fi

    docker_local "${IMAGE}" bash -lc '
      set -euo pipefail
      samtools faidx /local/ref/hg19.fa
      Rscript /work/bin/export_pace_regions.R --output /local/ref/pace.bed --flank 500
      Rscript /work/bin/export_pace_regions.R --output /local/ref/pace-screen.raw.bed --flank 1000
    '

    LC_ALL=C sort -k1,1 -k2,2n "${REF_DIR}/pace-screen.raw.bed" | \
      awk 'BEGIN{OFS="\t"} NR==1{c=$1;s=$2;e=$3;next} $1==c && $2<=e{if($3>e)e=$3;next} {print c,s,e;c=$1;s=$2;e=$3} END{if(NR>0)print c,s,e}' \
      > "${REF_DIR}/pace-screen.merged.bed"

    docker_local "${IMAGE}" bash -lc '
      set -euo pipefail
      : > /local/ref/pace-screen.fa
      while read -r chrom start end; do
        samtools faidx /local/ref/hg19.fa "${chrom}:$((start + 1))-${end}"
      done < /local/ref/pace-screen.merged.bed >> /local/ref/pace-screen.fa
      test -s /local/ref/pace-screen.fa
      samtools faidx /local/ref/pace-screen.fa
      bwameth.py index /local/ref/pace-screen.fa
      bwameth.py index /local/ref/hg19.fa
    '
    awk '{sum += $3-$2} END {print sum}' "${REF_DIR}/pace-screen.merged.bed" > "${REF_DIR}/pace-screen-reference-bases.txt"
    printf 'git_sha=%s\nimage=%s\ncreated=%s\n' "${GIT_SHA}" "${IMAGE}" "$(timestamp)" > "${ready}"
  } > "${ref_log}" 2>&1
  log "Referencja hg19 + PACE gotowa."
}

REF_PID=""
wait_for_reference() {
  if [[ -n "${REF_PID}" ]]; then
    log "Czekam na jednorazowe przygotowanie hg19/indeksow..."
    set +e
    wait "${REF_PID}"
    local status=$?
    set -e
    REF_PID=""
    if (( status != 0 )); then
      tail -n 80 "${LOG_DIR}/reference.log" >&2 || true
      die "Przygotowanie referencji nie powiodlo sie."
    fi
  fi
  [[ -s "${REF_DIR}/.hg19-pace-ready-v1" ]] || die "Brak markera gotowej referencji."
}

copy_compact_qc() {
  local sample="$1" out_dir="$2" sample_qc="$3" dest="${QC_DIR}/$1" suffix
  mkdir -p "${dest}"
  cp -f "${sample_qc}/"* "${dest}/" 2>/dev/null || true
  for suffix in pace-screen.source-stream.log pace-screen.fastp.json pace-screen.fastp.html \
    pace-screen.fastp.stderr.log pace-screen.bwameth.stderr.log pace-screen.samtools-view.stderr.log \
    pace-screen.samtools-fastq.stderr.log pace-screen.pipeline-status.txt source-stream.log \
    bwameth.stderr.log samtools-view.stderr.log pipeline-status.txt pace-targets.flagstat.txt; do
    cp -f "${out_dir}/${sample}.${suffix}" "${dest}/" 2>/dev/null || true
  done
}

process_sample() {
  local sample="$1" condition="$2" replicate="$3" r1_url="$4" r2_url="$5" run_accession="$6"
  local md5_r1="${7,,}" md5_r2="${8,,}"
  local sample_dir="${WORK_DIR}/${sample}" raw_dir="${WORK_DIR}/${sample}/raw"
  local out_dir="${WORK_DIR}/${sample}/out" sample_qc="${WORK_DIR}/${sample}/qc"
  local cov="${COV_DIR}/${sample}.bismark.cov.gz"
  local raw_r1="${WORK_DIR}/${sample}/raw/${sample}.R1.fastq.gz"
  local raw_r2="${WORK_DIR}/${sample}/raw/${sample}.R2.fastq.gz"
  local cand_r1="${WORK_DIR}/${sample}/out/${sample}.candidate.R1.fastq.gz"
  local cand_r2="${WORK_DIR}/${sample}/out/${sample}.candidate.R2.fastq.gz"
  local final_bam="${WORK_DIR}/${sample}/out/${sample}.mapq30.pace-targets.deduplicated.bam"
  local final_bai="${WORK_DIR}/${sample}/out/${sample}.mapq30.pace-targets.deduplicated.bam.bai"
  local start end p1 p2 s1 s2

  if [[ -s "${cov}" ]] && gzip -t "${cov}" >/dev/null 2>&1; then
    log "${sample}: coverage juz istnieje; SKIP."
    return 0
  fi

  check_disk
  mkdir -p "${raw_dir}" "${out_dir}" "${sample_qc}"

  if [[ ! -s "${final_bam}" || ! -s "${final_bai}" ]]; then
    if [[ ! -s "${cand_r1}" || ! -s "${cand_r2}" ]]; then
      log "${sample}: pobieram kompletne R1/R2 rownolegle z resume + published MD5."
      start="$(date +%s)"
      set +e
      download_with_resume "${sample} R1" "${r1_url}" "${raw_r1}" "${md5_r1}" \
        > >(tee -a "${sample_qc}/download.R1.log" "${RUN_LOG}") 2>&1 &
      p1=$!
      download_with_resume "${sample} R2" "${r2_url}" "${raw_r2}" "${md5_r2}" \
        > >(tee -a "${sample_qc}/download.R2.log" "${RUN_LOG}") 2>&1 &
      p2=$!
      wait "${p1}"; s1=$?
      wait "${p2}"; s2=$?
      set -e
      if (( s1 != 0 || s2 != 0 )); then
        log "${sample}: download nie skonczyl sie poprawnie (R1=${s1}, R2=${s2})."
        return 1
      fi
      end="$(date +%s)"
      echo $((end-start)) > "${sample_qc}/download-wall-seconds.txt"

      wait_for_reference
      log "${sample}: PACE prefilter na pelnych danych (screen +/-1000 bp)."
      start="$(date +%s)"
      if ! docker_local -e SAMPLE_ID="${sample}" -e MD5_R1="${md5_r1}" -e MD5_R2="${md5_r2}" -e THREADS="${THREADS}" \
        "${IMAGE}" bash -lc '
          set -euo pipefail
          cd "/local/work/${SAMPLE_ID}/out"
          bash /work/bin/screen_wgbs_pace_candidates.sh \
            "${SAMPLE_ID}" \
            "/local/work/${SAMPLE_ID}/raw/${SAMPLE_ID}.R1.fastq.gz" \
            "/local/work/${SAMPLE_ID}/raw/${SAMPLE_ID}.R2.fastq.gz" \
            "${MD5_R1}" "${MD5_R2}" /local/ref/pace-screen.fa "${THREADS}"
        ' 2>&1 | tee -a "${sample_qc}/prefilter.log" "${RUN_LOG}"; then
        log "${sample}: prefilter failed."
        return 1
      fi
      end="$(date +%s)"
      echo $((end-start)) > "${sample_qc}/prefilter-wall-seconds.txt"
      [[ -s "${cand_r1}" && -s "${cand_r2}" ]] || { log "${sample}: brak candidate FASTQ."; return 1; }
      gzip -t "${cand_r1}" && gzip -t "${cand_r2}" || return 1
      cp -f "${out_dir}/${sample}.candidate-pair-count.txt" "${sample_qc}/"
      cp -f "${out_dir}/${sample}.candidate-fastq-size.txt" "${sample_qc}/"
      rm -f "${raw_r1}" "${raw_r2}" "${raw_r1}.md5ok" "${raw_r2}.md5ok"
      log "${sample}: raw FASTQ usuniete po poprawnym prefilterze."
    else
      wait_for_reference
      log "${sample}: candidate FASTQ juz istnieja; wznawiam od full-hg19."
    fi

    log "${sample}: finalny alignment kandydatow do pelnego hg19, MAPQ >=30."
    start="$(date +%s)"
    if ! docker_local -e SAMPLE_ID="${sample}" -e THREADS="${THREADS}" -e WGBS_SKIP_FASTP=1 -e WGBS_SORT_MEM=512M \
      "${IMAGE}" bash -lc '
        set -euo pipefail
        cd "/local/work/${SAMPLE_ID}/out"
        bash /work/bin/stream_wgbs_bwameth.sh \
          "${SAMPLE_ID}" "${SAMPLE_ID}.candidate.R1.fastq.gz" "${SAMPLE_ID}.candidate.R2.fastq.gz" \
          "" "" /local/ref/hg19.fa /local/ref/pace.bed 30 "${THREADS}"
      ' 2>&1 | tee -a "${sample_qc}/full-hg19.log" "${RUN_LOG}"; then
      log "${sample}: full-hg19 alignment failed."
      return 1
    fi
    end="$(date +%s)"
    echo $((end-start)) > "${sample_qc}/full-hg19-candidate-wall-seconds.txt"
    [[ -s "${final_bam}" && -s "${final_bai}" ]] || { log "${sample}: brak finalnego BAM/BAI."; return 1; }
    rm -f "${cand_r1}" "${cand_r2}"
  else
    wait_for_reference
    log "${sample}: finalny BAM/BAI istnieje; wznawiam od methylation extraction."
  fi

  log "${sample}: MethylDackel + Bismark coverage."
  start="$(date +%s)"
  if ! docker_local -e SAMPLE_ID="${sample}" "${IMAGE}" bash -lc '
      set -euo pipefail
      cd "/local/work/${SAMPLE_ID}/out"
      micromamba run -n methyldackel MethylDackel extract -q 30 -p 5 --minDepth 1 -@ 4 \
        /local/ref/hg19.fa "${SAMPLE_ID}.mapq30.pace-targets.deduplicated.bam" -o "${SAMPLE_ID}.methyldackel"
      test -s "${SAMPLE_ID}.methyldackel_CpG.bedGraph"
      python3 /work/bin/methyldackel_to_bismark.py \
        --input "${SAMPLE_ID}.methyldackel_CpG.bedGraph" \
        --output "/local/results/cov/${SAMPLE_ID}.bismark.cov.gz"
      gzip -t "/local/results/cov/${SAMPLE_ID}.bismark.cov.gz"
      test "$(gzip -cd "/local/results/cov/${SAMPLE_ID}.bismark.cov.gz" | wc -l)" -gt 0
    ' 2>&1 | tee -a "${sample_qc}/methylation.log" "${RUN_LOG}"; then
    log "${sample}: methylation extraction failed."
    return 1
  fi
  end="$(date +%s)"
  echo $((end-start)) > "${sample_qc}/methylation-wall-seconds.txt"
  [[ -s "${cov}" ]] && gzip -t "${cov}" >/dev/null 2>&1 || return 1
  gzip -cd "${cov}" | wc -l > "${sample_qc}/coverage-row-count.txt"

  cat > "${sample_qc}/provenance.txt" <<EOF
sample_id=${sample}
condition=${condition}
replicate=${replicate}
run_accession=${run_accession}
source_scope=complete paired FASTQ archives, published MD5 verified
prefilter=fastp once; conservative DunedinPACE +/-1000bp screen; keep pair if either primary mate maps
final_reference=complete checksum-verified UCSC hg19
final_mapq_threshold=30
final_target_scope=official DunedinPACE probe windows +/-500bp
duplicate_removal=after complete-hg19 candidate alignment within biological sample
methylation_min_depth_raw=1
final_dunedinpace_min_depth=6
local_threads=${THREADS}
git_sha=${GIT_SHA}
image=${IMAGE}
EOF

  copy_compact_qc "${sample}" "${out_dir}" "${sample_qc}"
  rm -rf "${sample_dir}"
  log "${sample}: SUKCES. Coverage zachowane, ciezkie intermediates usuniete."
}

compute_pace_if_complete() {
  local expected="$1" actual
  actual="$(find "${COV_DIR}" -maxdepth 1 -type f -name '*.bismark.cov.gz' | wc -l | tr -d ' ')"
  log "Gotowe coverage: ${actual}/${expected}."
  (( actual == expected )) || { log "Nie licze jeszcze DunedinPACE: potrzebuje kompletu."; return 1; }

  log "Wszystkie coverage gotowe. Licze measured DunedinPACE."
  docker_local "${IMAGE}" bash -lc '
    set -euo pipefail
    Rscript /work/bin/compute_dunedinpace.R \
      --cov-dir /local/results/cov \
      --scores /local/results/dunedinpace_scores.csv \
      --qc /local/results/pace_probe_qc.csv \
      --model-metadata /local/results/pace_model_metadata.csv \
      --min-probe-fraction 0.80 --min-depth 6
    test -s /local/results/dunedinpace_scores.csv
    test -s /local/results/pace_probe_qc.csv
    test -s /local/results/pace_model_metadata.csv
  ' 2>&1 | tee -a "${RUN_LOG}"
  log "DunedinPACE GOTOWY: ${RESULTS_DIR}/dunedinpace_scores.csv"
}

trap 'log "Przerwano. Czesciowe dane zostaja; ponowne uruchomienie wznowi prace."' INT TERM

log "=== Bounty #1 local WGBS / MSYS2 ==="
log "Repo commit: ${GIT_SHA}"
log "Work root: ${LOCAL_ROOT}"
log "Threads: ${THREADS}; jedna probka naraz dla bezpieczenstwa 32 GB RAM."
check_disk
build_image

# Prepare hg19 in parallel with the first sample download; reuse it thereafter.
if [[ -s "${REF_DIR}/.hg19-pace-ready-v1" ]]; then
  log "Referencja jest juz przygotowana."
else
  prepare_reference &
  REF_PID=$!
  log "Przygotowanie referencji uruchomione w tle (PID ${REF_PID})."
fi

mapfile -t rows < <(
  awk -F',' 'NR>1 && $4=="WGBS" {gsub(/\r/,""); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$6,$7,$8,$11,$12}' "${SAMPLESHEET}"
)
selected=()
for row in "${rows[@]}"; do
  IFS=$'\t' read -r sample condition replicate r1_url r2_url run_accession md5_r1 md5_r2 <<< "${row}"
  [[ "${sample}" =~ ${SAMPLE_REGEX} ]] && selected+=("${row}")
done
(( ${#selected[@]} > 0 )) || die "SAMPLE_REGEX nie wybral zadnej probki."
log "Wybrano ${#selected[@]} WGBS samples."

failures=()
for row in "${selected[@]}"; do
  IFS=$'\t' read -r sample condition replicate r1_url r2_url run_accession md5_r1 md5_r2 <<< "${row}"
  log "----- START ${sample} -----"
  if process_sample "${sample}" "${condition}" "${replicate}" "${r1_url}" "${r2_url}" "${run_accession}" "${md5_r1}" "${md5_r2}"; then
    log "----- DONE ${sample} -----"
  else
    failures+=("${sample}")
    log "----- FAILED ${sample}; zachowuje stan do resume i ide dalej -----"
  fi
done

wait_for_reference
if (( ${#failures[@]} > 0 )); then
  log "Nieudane probki: ${failures[*]}"
  log "Pozostale probki zostaly przetworzone. Uruchom ten sam skrypt ponownie po diagnozie/retry."
  exit 1
fi

if (( ${#selected[@]} == ${#rows[@]} )); then
  compute_pace_if_complete "${#rows[@]}"
else
  log "Uruchomiono subset (${#selected[@]}/${#rows[@]}); nie agreguje DunedinPACE automatycznie."
fi
log "=== LOCAL WGBS FINISHED ==="
