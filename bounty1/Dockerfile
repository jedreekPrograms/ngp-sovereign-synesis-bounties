FROM mambaorg/micromamba:2.0.5

USER root

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.11 nextflow openjdk=21 fastp bowtie2 samtools macs3 bismark pysam scipy \
      r-base=4.3 r-optparse r-data.table r-remotes \
      bioconductor-minfi \
      bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19 \
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
