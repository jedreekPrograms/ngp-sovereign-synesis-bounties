# Bounty #1 — Reproducible ChIP-seq & DunedinPACE Pipeline

A reproducible Nextflow pipeline for the SIRT6 histone-acetylation / DunedinPACE bounty.

## Scientific data source

The default analysis is designed around open-access **HRA003336 / PRJCA012536** from the Institute of Zoology, Chinese Academy of Sciences. The project contains H3K9ac ChIP-seq, H3K56ac ChIP-seq and WGBS for WT and SIRT1–SIRT7-deficient human mesenchymal progenitor cells, including replicate samples.

The pipeline does **not** hard-code acceptance values. Correlations and p-values are computed only from processed sample-level outputs.

## Acceptance-relevant choices

- ChIP and WGBS BAM MAPQ filter: **30**
- MACS3 q-value threshold: **0.01**
- DunedinPACE: upstream `danbelsky/DunedinPACE` `PACEProjector()`
- DunedinPACE normalization: full `getRequiredProbes(backgroundList=TRUE)` probe set
- Pairing key: `condition + replicate`
- Histone occupancy: `log2(CPM + 1)` over union peaks, centered on the WT mean for each mark
- Pearson r and two-sided p-value: computed with SciPy after pairing ChIP occupancy with WGBS-derived DunedinPACE scores
- Paired-end and single-end FASTQ inputs are both supported

## Input samplesheet

Copy `resources/samplesheet.template.csv` to `resources/samplesheet.csv` and fill in the real FASTQ paths. Required columns:

```text
sample_id,condition,replicate,assay,mark,fastq_1,fastq_2
```

`fastq_2` may be blank for single-end data. ChIP `sample_id` values must contain the mark, for example `SIRT6_rep1_H3K9ac`; WGBS IDs should begin with the same `CONDITION_repN` key, for example `SIRT6_rep1_WGBS`.

## Run

Prepare hg19/GRCh37 Bowtie2 and Bismark indexes, then:

```bash
docker build -t bounty1-pace:1.0.0 .
docker run --rm \
  -v "$PWD:/work" \
  -v "/path/to/data:/data:ro" \
  -v "/path/to/refs:/refs:ro" \
  bounty1-pace:1.0.0 \
  nextflow run main.nf -resume \
  --samplesheet resources/samplesheet.csv \
  --bowtie2_index /refs/hg19/bowtie2/hg19 \
  --bismark_index /refs/hg19/bismark
```

## Outputs

The workflow creates:

- `results/qc/fastp/` — FASTQ QC
- `results/alignment/` — MAPQ-filtered BAMs
- `results/peaks/` — MACS3 peaks
- `results/methylation/` — Bismark CpG coverage
- `results/pace/dunedinpace_scores.csv` — sample-level DunedinPACE
- `results/chip/histone_occupancy.csv` — sample-level differential occupancy
- `results/correlations.json` — computed Pearson r/p values
- `results/manifest.json` — acceptance manifest generated from runtime outputs

## Integrity policy

This repository intentionally refuses to fabricate a DOI, correlation, p-value, FDR result, methylation score, or Docker checksum. `results/manifest.json` is generated from runtime outputs. A DOI and Docker-image checksum are added only when those artifacts actually exist.

## External data/model terms

The workflow does not vendor third-party raw data or DunedinPACE model assets. Users are responsible for complying with the access and licensing terms of HRA003336/NGDC and the upstream DunedinPACE package. The DunedinPACE project currently describes its algorithm as research-use-only for non-commercial users; obtain any license required for the intended use before a paid/commercial deployment.
