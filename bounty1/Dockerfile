FROM mambaorg/micromamba:2.0.5

USER root

# Exact versions are pinned to the set exercised by GitHub Actions on
# 2026-08-22. This prevents silent solver upgrades from changing scientific
# outputs between reruns.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.11.16 \
      nextflow=26.04.6 \
      openjdk=21.0.9 \
      fastp=1.3.6 \
      bowtie2=2.5.5 \
      bwa=0.7.19 \
      bwameth=0.2.9 \
      samtools=1.24 \
      macs3=3.0.4 \
      bismark=3.1.0 \
      pysam=0.24.0 \
      scipy=1.17.1 \
      matplotlib-base=3.10.9 \
      reportlab=5.0.0 \
      r-base=4.3.3 \
      r-optparse=1.7.5 \
      r-data.table=1.17.8 \
      r-remotes=2.5.0 \
      bioconductor-minfi=1.48.0 \
      bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19=0.6.0 \
    && micromamba clean --all --yes

# BSMAP is kept in an isolated environment because the Bioconda 2.90 build has
# legacy runtime dependencies. It is used only by the experimental, fail-closed
# WGBS prefilter benchmark; final production alignment remains bwa-meth against
# the complete hg19 reference.
RUN micromamba create -y -n bsmap -c conda-forge -c bioconda \
      bsmap=2.90 \
    && micromamba clean --all --yes

# MethylDackel 0.6.1 currently links against htslib <1.24, while the primary
# workflow intentionally pins samtools/htslib 1.24. Keep the extractor in an
# isolated environment rather than weakening the primary toolchain pin.
RUN micromamba create -y -n methyldackel -c conda-forge -c bioconda \
      methyldackel=0.6.1 \
    && micromamba clean --all --yes

# Pin the upstream research implementation for reproducibility. Docker build
# RUN steps do not execute through micromamba's runtime entrypoint, therefore
# R is invoked explicitly inside the base environment.
ARG DUNEDINPACE_REF=4b569983543e51d1022aecec9a25e694bb3a336a
RUN micromamba run -n base R -q -e \
    "remotes::install_github('danbelsky/DunedinPACE', ref='${DUNEDINPACE_REF}', upgrade='never')"

WORKDIR /work
COPY . /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["nextflow", "run", "main.nf"]
