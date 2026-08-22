# Bounty #1 acceptance-criteria audit

This document records acceptance requirements that were independently checked against the cited upstream implementation and the selected public dataset. It exists to keep the submission scientifically auditable and to distinguish measured outputs from values asserted by the bounty test harness.

## 1. DunedinPACE intercept mismatch

The bounty test `tests/test_bounty1_pace.py` requires:

```text
DunedinPACE intercept = 51.024577 +/- 0.001
```

The reproducible container in this submission pins the official research implementation:

```text
danbelsky/DunedinPACE
commit 4b569983543e51d1022aecec9a25e694bb3a336a
package version 0.99.0
```

At that pinned revision, `PACEProjector()` computes scores from `mPACE_Models$model_intercept` plus the weighted normalized methylation probes. The container CI directly loads that model object and reports:

```text
$DunedinPACE
[1] -1.949859
```

The official project documentation also describes a typical DunedinPACE result as approximately `1.0`, not approximately `51`.

Sources:

- https://github.com/danbelsky/DunedinPACE
- https://github.com/danbelsky/DunedinPACE/blob/main/R/PACEProjector.R
- Belsky et al. (2022), eLife 11:e73420, https://doi.org/10.7554/eLife.73420

A search of the official implementation/documentation and publication material did not identify `51.024577` as a DunedinPACE model intercept. Consequently, this repository **does not hard-code 51.024577**. The runtime manifest reads the intercept from the pinned upstream model and therefore currently records `-1.949859`.

If the bounty issuer intended a different historical model, preprocessing transformation, or regression intercept, that model must be identified before the numerical criterion can be reproduced honestly.

## 2. WGBS versus methylation-array input

The bounty explicitly permits paired WGBS or 450K/EPIC methylation data. The official DunedinPACE software itself is documented for Illumina 450K/EPIC/EPICv2 beta matrices.

For WGBS, this pipeline therefore performs a transparent coordinate projection:

1. obtain the full DunedinPACE background-probe list with `getRequiredProbes(backgroundList=TRUE)`;
2. use the EPIC hg19 annotation to map those `cg...` probes to genomic CpG coordinates;
3. extract WGBS beta values only at those coordinates;
4. require depth >= 5 reads at a CpG;
5. require >= 80% of required probes for every sample;
6. pass the resulting beta matrix to the unmodified upstream `PACEProjector()`.

This is a WGBS-to-array-probe projection, not a claim that the upstream package natively consumes Bismark coverage files.

## 3. Public dataset and biological replication

The production metadata comes from the official HRA003336 Excel export and is tracked in `resources/samplesheet.csv` and `resources/chip_inputs.csv`.

Primary analysis matrix:

- 16 H3K9ac ChIP-seq libraries;
- 16 H3K56ac ChIP-seq libraries;
- 16 WGBS libraries;
- 16 fully matched `condition + replicate` triples;
- WT plus SIRT1-SIRT7;
- two biological replicates per condition.

Controls:

- eight condition-matched ChIP INPUT libraries, one for each of WT and SIRT1-SIRT7.

The objective text says `across >= 3 biological replicates`. The selected dataset has 16 paired biological sample records overall but only two replicates within each condition. The submission reports that structure explicitly rather than relabeling samples to manufacture a third replicate.

## 4. Correlation criterion

The required Pearson thresholds (`r > 0.92`, `p < 0.01`) are treated as hypotheses/acceptance targets, not constants. `bin/correlate.py` calculates them from sample-level ChIP occupancy and WGBS-derived DunedinPACE outputs. If the real data do not meet the threshold, the pipeline reports the measured result and fails the criterion rather than substituting a synthetic value.

## 5. Other acceptance settings

The following settings are direct, reproducible pipeline parameters and already satisfy the numerical requirements by construction:

- BAM MAPQ threshold: `30` (required >=30)
- MACS3 q-value threshold: `0.01` (required <0.05)
- condition-matched INPUT control used for MACS3 peak calling
- hg19/GRCh37 coordinates throughout the ChIP/WGBS branches and EPIC-probe projection

## 6. DOI and Docker checksum

The pipeline will not invent either artifact.

- `data_deposit_doi` remains empty until a real archival deposit exists.
- `docker_image_md5` is emitted only for an actually exported Docker image tarball.

These are final-release tasks after production processing and report generation.

## 7. Measured full-source technical benchmark evidence

A complete-source SIRT6 CUT&RUN technical benchmark completed successfully in GitHub Actions run `32595547809` on 2026-08-22. This benchmark intentionally aligns against UCSC hg19 chromosome 22 only, so it is scalability/integrity evidence and **not** a final whole-genome occupancy result.

Measured inputs and QC:

- treatment: `WT_rep2_SIRT6_CUTRUN` (`HRR1354864`), 16,813,112 paired FASTQ records before filtering;
- matched control: `WT_rep2_IgG_CUTRUN` (`HRR1354852`), 12,498,528 paired FASTQ records before filtering;
- all four compressed-source MD5 values matched the manifest before the samples were accepted;
- treatment retained 176,024 MAPQ >= 30 reads after alignment, duplicate removal, and filtering on chr22;
- control retained 159,428 MAPQ >= 30 reads after the same processing;
- MACS3 was run with `q <= 0.05` and produced 882 chr22 peaks;
- measured wall time was 1,844 s for treatment and 1,519 s for control;
- final benchmark working set was approximately 130 MiB after intermediate cleanup.

The uploaded evidence artifact is `bounty1-full-source-cutrun-chr22` with artifact ZIP SHA-256:

```text
f6c8a09f65324534f9c86a002120e44ccff60a59a9f2c38d84e16f06ad16bc32
```

The artifact contains fastp JSON/HTML reports, Bowtie2 logs, MAPQ30 flagstats, MACS3 logs, the narrowPeak file, reference checksum, wall-time measurements, disk snapshots, and benchmark-scope documentation.

## 8. Measured full-source H3K9ac ChIP technical benchmark evidence

A complete-source SIRT7 H3K9ac ChIP technical benchmark also completed successfully in GitHub Actions run `32595547809` on 2026-08-22. Like the CUT&RUN benchmark above, this intentionally aligns against UCSC hg19 chromosome 22 only and is therefore scalability/integrity evidence rather than a final whole-genome occupancy result.

Measured inputs and QC:

- treatment: `SIRT7_rep1_H3K9ac` (`HRR1202646`), 38,740,454 paired FASTQ records before filtering;
- matched INPUT control: `SIRT7_INPUT` (`HRR1202672`), 30,746,934 paired FASTQ records before filtering;
- all four compressed-source MD5 values matched the manifest before the samples were accepted;
- treatment retained 604,496 MAPQ >= 30 reads after alignment, duplicate removal, and filtering on chr22;
- control retained 623,308 MAPQ >= 30 reads after the same processing;
- MACS3 was run with `q <= 0.01` and produced 909 chr22 peaks;
- measured wall time was 3,620 s for treatment and 4,457 s for control;
- final benchmark working set was approximately 242 MiB after intermediate cleanup.

The uploaded evidence artifact is `bounty1-full-source-chip-chr22` with artifact ZIP SHA-256:

```text
0b6da3b6867be91dfc36503dd5a918ae011dbe2b064cf4c34d0a5c8c66815c1b
```

The artifact contains fastp JSON/HTML reports, Bowtie2 logs, MAPQ30 flagstats, MACS3 logs, the narrowPeak file, reference checksum, wall-time measurements, disk snapshots, and benchmark-scope documentation.
