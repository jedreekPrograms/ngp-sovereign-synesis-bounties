# Bounty #1 remote-compute escape plan

The GitHub-hosted WGBS cohort remains useful for benchmarking and checkpoints, but production should not depend on a single <=6 h hosted job for each 40-55 GiB WGBS sample.

## Priority order

1. **Obtain authentic processed WGBS methylation output if available.** The bounty explicitly requires raw FASTQ for ChIP-seq, while its methylation objective allows paired WGBS or 450K/EPIC methylation data. The HRA003336 paper states that the sequencing data are public and names Dr. Jing Qu (qujing@ioz.ac.cn) as the lead contact for additional information/resources. Request per-sample, pre-replicate-merge hg19 CpG methylation output (preferably counts/depth plus methylation ratio) for all 16 WGBS samples. If supplied, preserve provenance and checksums and do not invent missing values. This could remove most of the ~740 GiB WGBS raw-compute burden.

2. **Free academic HPC: WCSS / PCSS.** Run the existing Nextflow workflow on SLURM using `config/hpc-slurm.config`. Stage source FASTQs to project scratch with resumable transfers and published MD5 verification, then let Nextflow resume from its work directory. Do not use GitHub artifacts for raw FASTQ storage.

3. **Managed cloud: Seqera Compute.** The project is already Nextflow DSL2. A managed environment close to the CNCB source (for example AWS `ap-east-1`, Hong Kong, when available to the account) avoids the GitHub-hosted 6 h job ceiling. Keep raw data streaming/ephemeral and persist only validated reduced outputs. Benchmark one WGBS sample first and set a hard spend ceiling before launching the full cohort.

4. **Ephemeral self-hosted GitHub runner.** Attach a remote Linux VM with large local SSD to this repository and keep GitHub as the orchestration UI. Use a dedicated runner label and manual-only production workflow. Destroy/deregister the runner after the production run. Never expose a persistent self-hosted runner to untrusted pull-request code in this public repository.

5. **GitHub-hosted fallback.** Continue checkpointed v3 and, if needed, replace the current monolithic bisulfite prefilter with a much faster conservative candidate screen. Any new screening method must be measured against a direct authentic-data alignment subset for recall before being used as final production evidence.

## Data-integrity rules

- Published source MD5 must be checked for every complete source FASTQ used as evidence.
- Intermediate checkpoints get SHA-256 plus provenance metadata.
- Final WGBS candidate alignment remains against complete checksum-verified hg19 with MAPQ >= 30.
- Do not split gzip files at arbitrary compressed byte offsets unless an indexing/container strategy proves independent decompression and preserved R1/R2 record pairing.
- Do not merge biological replicates before producing the per-sample DunedinPACE inputs required for the correlation analysis.
- Never substitute fixture/stub values for measured biological outputs.

## HPC invocation skeleton

```bash
export NXF_WORK=/path/to/project/scratch/bounty1-work
nextflow -c bounty1/nextflow.config \
  -c bounty1/config/hpc-slurm.config \
  run bounty1/main.nf \
  -resume
```

Site-specific SLURM account/partition directives should live in a private local config rather than being guessed in the public repository.
