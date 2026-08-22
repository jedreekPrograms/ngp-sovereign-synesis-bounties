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
      samtools=1.24 \
      macs3=3.0.4 \
      bismark=3.1.0 \
      pysam=0.24.0 \
      scipy=1.17.1 \
      r-base=4.3.3 \
      r-optparse=1.7.5 \
      r-data.table=1.17.8 \
      r-remotes=2.5.0 \
      bioconductor-minfi=1.48.0 \
      bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19=0.6.0 \
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
