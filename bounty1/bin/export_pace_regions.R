#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DunedinPACE)
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option('--output', type='character', default='dunedinpace_probe_windows.bed'),
  make_option('--flank', type='integer', default=500)
)))

if (opts$flank < 0) stop('--flank must be >= 0')

required <- unique(unlist(getRequiredProbes(backgroundList=TRUE), use.names=FALSE))
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- ann[intersect(required, rownames(ann)), , drop=FALSE]
if (nrow(ann) == 0) stop('No DunedinPACE probes found in EPIC hg19 annotation')

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl('^chr', x), x, paste0('chr', x))
}

# EPIC annotation positions are 1-based. BED is 0-based, half-open. Retaining a
# flank around each probe makes it overwhelmingly likely that both mates of a
# fragment are preserved before duplicate marking while still reducing the
# WGBS alignment footprint to a small fraction of the genome.
pos0 <- pmax(as.integer(ann$pos) - 1L, 0L)
windows <- data.table(
  chr = normalize_chr(ann$chr),
  start = pmax(pos0 - opts$flank, 0L),
  end = pos0 + opts$flank + 1L,
  probe = rownames(ann)
)
windows <- unique(windows[!is.na(chr) & !is.na(start) & !is.na(end)])
setorder(windows, chr, start, end, probe)

if (nrow(windows) < 1000) {
  stop(sprintf('Unexpectedly small DunedinPACE region set: %d probes', nrow(windows)))
}

fwrite(windows, opts$output, sep='\t', col.names=FALSE)
message(sprintf('Wrote %d DunedinPACE probe windows to %s', nrow(windows), opts$output))
