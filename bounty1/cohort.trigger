Launch the full 16-sample HRA003336 WGBS cohort using the conservative PACE-candidate prefilter workflow.

Triggered after the direct full-source SIRT7_rep1 benchmark demonstrated that unrestricted full-genome bwa-meth exceeds the standard GitHub-hosted runner time budget. Final candidate alignments still use complete checksum-verified hg19 and MAPQ >= 30; the prefilter only reduces the reads sent to that expensive step.
