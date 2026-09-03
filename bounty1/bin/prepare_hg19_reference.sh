#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-refs/hg19}"
UCSC_BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips"
FASTA_GZ="hg19.fa.gz"

mkdir -p "${OUTDIR}/download" "${OUTDIR}/bowtie2" "${OUTDIR}/bismark"

if [[ ! -s "${OUTDIR}/download/${FASTA_GZ}" ]]; then
  curl --fail --location --retry 5 --retry-delay 5 \
    "${UCSC_BASE}/${FASTA_GZ}" \
    --output "${OUTDIR}/download/${FASTA_GZ}"
fi

# Verify against the checksum published by UCSC instead of trusting transport.
curl --fail --location --retry 5 \
  "${UCSC_BASE}/md5sum.txt" \
  --output "${OUTDIR}/download/md5sum.txt"
(
  cd "${OUTDIR}/download"
  grep -E "[[:space:]]${FASTA_GZ}$" md5sum.txt > hg19.fa.gz.md5
  test -s hg19.fa.gz.md5
  md5sum --check hg19.fa.gz.md5
)

if [[ ! -s "${OUTDIR}/hg19.fa" ]]; then
  gzip -dc "${OUTDIR}/download/${FASTA_GZ}" > "${OUTDIR}/hg19.fa"
fi

# Regular Bowtie2 reference for ChIP-seq.
if [[ ! -s "${OUTDIR}/bowtie2/hg19.1.bt2" && ! -s "${OUTDIR}/bowtie2/hg19.1.bt2l" ]]; then
  bowtie2-build --threads "${BOWTIE2_BUILD_THREADS:-4}" \
    "${OUTDIR}/hg19.fa" "${OUTDIR}/bowtie2/hg19"
fi

# Bismark requires the FASTA to live in the genome folder it prepares.
cp -f "${OUTDIR}/hg19.fa" "${OUTDIR}/bismark/hg19.fa"
if [[ ! -d "${OUTDIR}/bismark/Bisulfite_Genome" ]]; then
  bismark_genome_preparation --bowtie2 --parallel "${BISMARK_BUILD_THREADS:-4}" \
    "${OUTDIR}/bismark"
fi

{
  sha256sum "${OUTDIR}/download/${FASTA_GZ}"
  sha256sum "${OUTDIR}/hg19.fa"
  find "${OUTDIR}/bowtie2" "${OUTDIR}/bismark/Bisulfite_Genome" \
    -type f -print0 | sort -z | xargs -0 sha256sum
} > "${OUTDIR}/reference_checksums.sha256"

cat <<EOF
Reference prepared successfully.
Bowtie2 prefix: ${OUTDIR}/bowtie2/hg19
Bismark genome: ${OUTDIR}/bismark
Checksums:      ${OUTDIR}/reference_checksums.sha256
EOF
