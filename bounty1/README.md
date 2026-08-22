# Bounty #1 — Reproducible ChIP-seq & DunedinPACE Pipeline

A reproducible Nextflow pipeline for the SIRT6 histone-acetylation / DunedinPACE bounty. The implementation is deliberately evidence-first: production correlations, p-values, DOI values and container checksums are generated from real artifacts rather than embedded acceptance constants.

## Scientific data sources

### Primary paired ChIP-seq + WGBS cohort

The default analysis uses open-access **HRA003336 / PRJCA012536** from the Institute of Zoology, Chinese Academy of Sciences. The study contains H3K9ac ChIP-seq, H3K56ac ChIP-seq and WGBS for WT and SIRT1–SIRT7-deficient human mesenchymal progenitor cells (hMPCs).

`resources/samplesheet.csv` was derived from the official `HRA003336.xlsx` metadata export and contains **48 paired-end assay libraries**:

- 16 H3K9ac ChIP-seq libraries,
- 16 H3K56ac ChIP-seq libraries,
- 16 WGBS libraries,
- WT + SIRT1–SIRT7,
- two biological replicates per condition.

`resources/chip_inputs.csv` adds **8 condition-matched ChIP INPUT libraries**.

### Independent SIRT6 CUT&RUN sensitivity dataset

`resources/sirt6_cutrun.csv` contains **8 paired-end HRA005392 libraries** from the same publication:

- WT and SIRT6-KO,
- two replicates per condition,
- SIRT6 antibody and matched IgG controls.

The SIRT6 CUT&RUN branch is a **prespecified secondary mechanistic analysis**. It does not replace the primary global histone-occupancy endpoint after results are observed.

### Archive inventory

The production manifests therefore reference **64 paired-end libraries / 128 FASTQ objects**. GitHub Actions probes every URL without downloading the payload. The current archive inventory reports:

- reachable FASTQ objects: **128 / 128**,
- HTTP/size failures: **0**,
- combined compressed payload: **1005.79 GiB**.

The CI planner also ranks complete `H3K9ac + H3K56ac + WGBS` observations by transfer size. The three smallest complete observations alone total approximately **155 GiB**, which is why a bounded real-data smoke test is kept separate from the full production run.

Official sources:

- https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003336
- https://ngdc.cncb.ac.cn/gsa-human/browse/HRA005392
- https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA012536
- Bi et al., *Developmental Cell* (2024), DOI: `10.1016/j.devcel.2024.02.008`

The paper's STAR Methods process these data against **hg19**, so the workflow deliberately uses hg19/GRCh37 for ChIP-seq, CUT&RUN, WGBS and EPIC probe-coordinate mapping.

## Prespecified analysis plan

### Primary bounty endpoint

For each histone mark separately:

1. call MACS3 peaks from the real ChIP-seq libraries using the condition-matched INPUT,
2. construct one fixed union peak universe for the mark,
3. count unique DNA fragments for every sample over exactly that same universe,
4. convert to `log2(fragment CPM + 1)`,
5. subtract the WT mean for that mark,
6. pair occupancy with WGBS-derived DunedinPACE by `condition + replicate`,
7. compute Pearson `r` and two-sided `p` with SciPy.

The peak universe is created without access to DunedinPACE values. `differential_occupancy` in `histone_occupancy.csv` is always this primary global endpoint.

### Secondary SIRT6-locus endpoint

SIRT6 CUT&RUN is processed independently. Reproducible WT SIRT6 peaks are required to overlap between the two WT replicates by at least 50% reciprocally; reproducible SIRT6-KO antibody peaks are then removed. H3K9ac/H3K56ac fragment occupancy is also quantified over these independent SIRT6-specific loci and reported as a secondary sensitivity analysis.

Both analyses are emitted. The secondary analysis is never silently substituted for the primary endpoint.

## Acceptance-relevant processing choices

- ChIP-seq and CUT&RUN MAPQ filter: **>= 30**
- WGBS MAPQ filter: **>= 30 for both mates of a retained pair**
- H3K9ac/H3K56ac MACS3 q-value: **0.01**
- SIRT6 CUT&RUN MACS3 q-value: **0.05**, matching the source-study methods
- PCR duplicate removal before final MAPQ filtering for ChIP-seq and CUT&RUN
- Bismark deduplication for WGBS
- pair-aware WGBS MAPQ filtering so R1/R2 cannot be orphaned
- paired WGBS BAM remains query-name sorted for Bismark methylation extraction
- reference: **hg19/GRCh37**
- WGBS CpG depth threshold: **>= 6 reads** (source paper specifies depth >5)
- DunedinPACE implementation pinned to upstream commit `4b569983543e51d1022aecec9a25e694bb3a336a`
- normalization uses `getRequiredProbes(backgroundList=TRUE)`, including the full background set
- DunedinPACE model intercept is read from the installed upstream model at runtime

The repository does **not** contain the bounty's expected intercept constant in production or unit-test fixtures. The mismatch between the bounty acceptance test and the pinned upstream DunedinPACE model is documented in `docs/acceptance_audit.md` and has been raised publicly on upstream Issue #1.

## Input provenance and integrity

Every production row records:

```text
sample_id,condition,replicate,assay,mark,fastq_1,fastq_2,
run_accession,experiment_accession,gsa_sample_accession,read1_md5,read2_md5
```

`bin/validate_samplesheet.py` validates the assay matrix, accession syntax, run uniqueness, FASTQ naming and MD5 syntax.

The HRA metadata exports contain a duplicated R1 URL in some displayed R2 URL fields. R2 URLs in this repository are reconstructed from the authoritative run accession + R2 filename while retaining the independent archive-provided R2 MD5.

## Streaming full raw data

Because the public payload is roughly 1 TiB compressed, production execution does not first stage the entire raw archive on local disk. `bin/stream_fastp.sh` streams compressed bytes from each public URL while calculating the archive MD5, decompresses the same byte stream, and feeds it to `fastp`. Only trimmed FASTQs and downstream products enter the Nextflow work tree.

`stream_io` limits concurrent archive streams so a local or cloud runner does not create uncontrolled network and disk pressure.

Local pre-downloaded FASTQ files are supported by the same helper.

## WGBS pair integrity

A naïve per-read `samtools view -q 30` can retain one mate and drop the other. That is unsafe for paired-end Bismark methylation extraction.

The production WGBS branch therefore:

1. aligns with Bismark,
2. deduplicates with `deduplicate_bismark`,
3. query-name sorts the deduplicated BAM,
4. runs `bin/filter_paired_bam.py`,
5. retains a fragment only when exactly one primary R1 and one primary R2 are present and **both** satisfy MAPQ >= 30,
6. preserves query-name order for `bismark_methylation_extractor`.

A JSON filter report is emitted per WGBS sample.

## DunedinPACE from WGBS

`bin/compute_dunedinpace.R` maps Bismark CpG coverage coordinates onto the hg19 Illumina EPIC annotation, constructs the beta matrix required by the upstream package and calls the pinned upstream `PACEProjector()`.

The full normalization background returned by `getRequiredProbes(backgroundList=TRUE)` is used. The model intercept recorded in the final manifest is read from `mPACE_Models$model_intercept` inside the installed package rather than typed into this repository.

**Scientific limitation:** DunedinPACE was developed primarily for blood DNA methylation and the upstream package documents Illumina 450K/EPIC beta-matrix input. Applying it to hMPC WGBS projected onto those CpG coordinates is an off-domain research use and will be stated explicitly in the final report.

## Reference preparation

`bin/prepare_hg19_reference.sh` prepares the hg19 Bowtie2 and Bismark reference assets and records checksums.

Example inside the analysis image:

```bash
bash bin/prepare_hg19_reference.sh /refs/hg19
```

## Production run

```bash
docker build -t bounty1-pace:1.0.0 .
docker run --rm \
  -v "$PWD:/work" \
  -v "/path/to/refs:/refs:ro" \
  bounty1-pace:1.0.0 \
  nextflow run main.nf -resume \
  --samplesheet resources/samplesheet.csv \
  --controls_sheet resources/chip_inputs.csv \
  --cutrun_sheet resources/sirt6_cutrun.csv \
  --bowtie2_index /refs/hg19/bowtie2/hg19 \
  --bismark_index /refs/hg19/bismark
```

## Bounded real-data smoke test

`.github/workflows/bounty1-real-smoke.yml` performs a technical smoke test on **real public reads**, not generated FASTQ:

- HRA003336 ChIP-seq + INPUT,
- HRA005392 SIRT6 CUT&RUN + IgG,
- HRA003336 WGBS,
- UCSC hg19 chromosome 22,
- first 75,000 paired records per selected library.

It exercises real FASTQ parsing, fastp, Bowtie2, deduplication, MAPQ filtering, MACS3 invocation, Bismark alignment, pair-aware WGBS filtering and methylation extraction. The workflow labels its outputs as **partial technical evidence only**. They are never used as final biological correlations.

## Main outputs

- `results/qc/fastp/` — FASTQ QC JSON/HTML
- `results/alignment/chip/` — deduplicated MAPQ-filtered ChIP BAMs
- `results/alignment/cutrun/` — deduplicated MAPQ-filtered CUT&RUN BAMs
- `results/alignment/wgbs/` — pair-preserving MAPQ-filtered WGBS BAMs + filter stats
- `results/peaks/` — histone MACS3 peaks
- `results/cutrun/` — CUT&RUN peaks and SIRT6-specific locus set
- `results/methylation/` — Bismark CpG coverage
- `results/pace/dunedinpace_scores.csv` — sample-level DunedinPACE
- `results/pace/pace_probe_qc.csv` — required-probe coverage per WGBS sample
- `results/chip/histone_occupancy.csv` — primary global + secondary SIRT6-locus occupancy
- `results/correlations.json` — primary and secondary Pearson analyses
- `results/manifest.json` — runtime-generated acceptance/provenance manifest

## Reproducibility and integrity policy

The Dockerfile pins the analysis stack and the exact DunedinPACE Git commit. GitHub Actions builds the real image, validates the production metadata, loads the installed DunedinPACE model, exercises the streaming helper, performs a full Nextflow `-stub-run`, checks every public archive URL and records a reproducible data-size inventory.

This repository intentionally refuses to fabricate:

- a correlation,
- a p-value,
- a DunedinPACE model intercept,
- a DOI,
- a Docker checksum,
- a biological result.

A DOI and final container checksum are added only after the corresponding real artifacts exist.

## External software/data terms

The workflow does not vendor the public human sequencing payloads or DunedinPACE model files. Users must comply with the source archives and third-party software terms.

The current DunedinPACE repository contains a GPL-3.0 `LICENSE` file while its README separately states that the algorithm/site is for research users and directs commercial users to the exclusive licensee. This bounty implementation is being developed as a public research/reproducibility analysis, not as a commercial diagnostic product. Any separate commercial deployment should resolve the upstream licensing notice directly with the rights holder/licensee.
