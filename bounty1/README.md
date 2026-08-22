# Bounty #1 — Reproducible ChIP-seq & DunedinPACE Pipeline

A reproducible Nextflow pipeline for the SIRT6 histone-acetylation / DunedinPACE bounty.

## Scientific data source

The default analysis uses open-access **HRA003336 / PRJCA012536** from the Institute of Zoology, Chinese Academy of Sciences. The project contains H3K9ac ChIP-seq, H3K56ac ChIP-seq and WGBS for WT and SIRT1–SIRT7-deficient human mesenchymal progenitor cells.

The exact run mapping in `resources/samplesheet.csv` was derived from the official `HRA003336.xlsx` metadata export. It contains **48 paired-end libraries**: 16 H3K9ac, 16 H3K56ac and 16 WGBS libraries, covering WT + SIRT1–SIRT7 and two replicates for each condition. Each row records the GSA sample (`HRS`), experiment (`HRX`), run (`HRR`), both FASTQ URLs and both archive-provided MD5 checksums. `resources/chip_inputs.csv` adds eight condition-matched INPUT libraries.

A CI availability probe checks all **112 FASTQ objects** (56 paired-end libraries) without downloading their payloads. At the time of validation all 112 URLs were reachable and exposed a combined compressed size of **988.15 GiB**.

Official sources:

- https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003336
- https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA012536
- Bi et al., Developmental Cell (2024), DOI: 10.1016/j.devcel.2024.02.008

The paper's STAR Methods process the ChIP-seq and WGBS data against **hg19**, so this workflow deliberately uses hg19/GRCh37 for both assay branches and for the EPIC probe-coordinate mapping used to construct the DunedinPACE beta matrix.

The pipeline does **not** hard-code acceptance values. Correlations and p-values are computed only from processed sample-level outputs. See `docs/acceptance_audit.md` for the independently identified DunedinPACE-intercept discrepancy in the bounty specification.

## Acceptance-relevant choices

- ChIP and WGBS BAM MAPQ filter: **30**
- MACS3 q-value threshold: **0.01**
- Condition-matched INPUT control for ChIP peak calling
- Reference genome: **hg19/GRCh37**, matching the source study
- DunedinPACE: pinned upstream `danbelsky/DunedinPACE` `PACEProjector()`
- DunedinPACE normalization: full `getRequiredProbes(backgroundList=TRUE)` probe set
- WGBS CpG depth threshold: **>= 6 reads** (source paper specifies depth >5)
- Pairing key: `condition + replicate`
- Histone occupancy: `log2(CPM + 1)` over union peaks, centered on the WT mean for each mark
- Pearson r and two-sided p-value: computed with SciPy after pairing ChIP occupancy with WGBS-derived DunedinPACE scores
- Paired-end and single-end FASTQ inputs are both supported

## Input samplesheet

A production samplesheet is already provided at `resources/samplesheet.csv`. Required execution columns are:

```text
sample_id,condition,replicate,assay,mark,fastq_1,fastq_2
```

The production sheet also carries provenance/integrity columns:

```text
run_accession,experiment_accession,gsa_sample_accession,read1_md5,read2_md5
```

`bin/validate_samplesheet.py` verifies the assay matrix, unique HRR runs, accession syntax, FASTQ naming and MD5 syntax before analysis. The GSA metadata export contains a duplicated R1 URL in its displayed "DownLoad Read file2" field; this repository therefore derives each R2 URL from the official HRR directory plus the independent `Read filename 2` value and preserves the archive-provided R2 MD5.

### Streaming raw data

Because the selected archive payload is nearly 1 TiB compressed, production execution does **not** ask Nextflow to stage every raw FASTQ onto local disk. `bin/stream_fastp.sh` streams each remote compressed FASTQ directly into `fastp` through FIFOs while simultaneously hashing the exact compressed byte stream and comparing it with the archive-provided MD5. Only trimmed FASTQs and downstream products enter the workflow work directory. The default `stream_io` process label limits concurrent archive streams to avoid uncontrolled network and disk pressure.

Local FASTQ paths are supported by the same streaming helper, so a pre-downloaded dataset can be substituted without changing the workflow logic.

## Reference preparation

`bin/prepare_hg19_reference.sh` downloads UCSC hg19, verifies the UCSC-published MD5, constructs the regular Bowtie2 index and Bismark bisulfite genome, and records SHA-256 checksums for the generated reference assets.

Example inside the analysis image:

```bash
bash bin/prepare_hg19_reference.sh /refs/hg19
```

## Run

Build the pinned image and run Nextflow:

```bash
docker build -t bounty1-pace:1.0.0 .
docker run --rm \
  -v "$PWD:/work" \
  -v "/path/to/refs:/refs:ro" \
  bounty1-pace:1.0.0 \
  nextflow run main.nf -resume \
  --samplesheet resources/samplesheet.csv \
  --controls_sheet resources/chip_inputs.csv \
  --bowtie2_index /refs/hg19/bowtie2/hg19 \
  --bismark_index /refs/hg19/bismark
```

## Outputs

The workflow creates:

- `results/qc/fastp/` — compact FASTQ QC JSON/HTML (trimmed FASTQs stay in the work tree rather than being duplicated into results)
- `results/alignment/` — MAPQ-filtered BAMs
- `results/peaks/` — MACS3 peaks
- `results/methylation/` — Bismark CpG coverage
- `results/pace/dunedinpace_scores.csv` — sample-level DunedinPACE
- `results/pace/pace_probe_qc.csv` — required-probe coverage per WGBS sample
- `results/chip/histone_occupancy.csv` — sample-level differential occupancy
- `results/correlations.json` — computed Pearson r/p values
- `results/manifest.json` — acceptance manifest generated from runtime outputs

## Reproducibility

The Dockerfile pins the exact analysis stack (Nextflow, Java, fastp, Bowtie2, samtools, MACS3, Bismark, Python/R packages) and pins DunedinPACE to a specific upstream Git commit. GitHub Actions builds the real image, validates production metadata, loads the actual DunedinPACE model, performs a complete Nextflow `-stub-run`, and checks all archive URLs.

## Integrity policy

This repository intentionally refuses to fabricate a DOI, correlation, p-value, FDR result, methylation score, or Docker checksum. `results/manifest.json` is generated from runtime outputs. A DOI and Docker-image checksum are added only when those artifacts actually exist.

## External data/model terms

The workflow does not vendor third-party raw data or DunedinPACE model assets. Users are responsible for complying with the access and licensing terms of HRA003336/NGDC and the upstream DunedinPACE package. The DunedinPACE project currently describes its algorithm as research-use-only for non-commercial users; obtain any license required for the intended use before a paid/commercial deployment.
