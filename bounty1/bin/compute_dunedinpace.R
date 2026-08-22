#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DunedinPACE)
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option('--cov-dir', type='character', default='.'),
  make_option('--scores', type='character', default='dunedinpace_scores.csv'),
  make_option('--qc', type='character', default='pace_probe_qc.csv'),
  make_option('--model-metadata', type='character', default='pace_model_metadata.csv'),
  make_option('--min-probe-fraction', type='double', default=0.80)
)))

required <- unique(unlist(getRequiredProbes(backgroundList=TRUE), use.names=FALSE))
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- ann[intersect(required, rownames(ann)), , drop=FALSE]
if (nrow(ann) == 0) stop('No DunedinPACE probes found in EPIC hg19 annotation')

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl('^chr', x), x, paste0('chr', x))
}

probe_table <- data.table(
  probe = rownames(ann),
  chr = normalize_chr(ann$chr),
  pos = as.integer(ann$pos)
)
setkey(probe_table, chr, pos)

cov_files <- list.files(opts$`cov-dir`, pattern='\\.bismark\\.cov\\.gz$', full.names=TRUE)
if (length(cov_files) < 3) stop('Need at least 3 WGBS samples')

beta <- matrix(
  NA_real_, nrow=length(required), ncol=length(cov_files),
  dimnames=list(required, sub('\\.bismark\\.cov\\.gz$', '', basename(cov_files)))
)
qc <- list()

for (f in cov_files) {
  sample_id <- sub('\\.bismark\\.cov\\.gz$', '', basename(f))
  x <- fread(f, col.names=c('chr','start','end','pct','meth','unmeth'))
  x[, chr := normalize_chr(chr)]
  x[, pos := as.integer(start)]
  x[, depth := meth + unmeth]
  x[, beta := meth / pmax(depth, 1)]
  x <- x[depth >= 5, .(chr, pos, beta)]
  setkey(x, chr, pos)

  joined <- x[probe_table, on=.(chr,pos), nomatch=0]
  joined <- joined[!is.na(probe)]
  vals <- joined[, .(beta=mean(beta, na.rm=TRUE)), by=probe]
  beta[vals$probe, sample_id] <- vals$beta

  matched <- sum(!is.na(beta[, sample_id]))
  coverage <- matched / length(required)
  qc[[sample_id]] <- data.table(
    sample_id=sample_id,
    matched_probes=matched,
    required_probes=length(required),
    fraction=coverage
  )
  if (coverage < opts$`min-probe-fraction`) {
    stop(sprintf('%s probe coverage %.3f below threshold %.3f', sample_id, coverage, opts$`min-probe-fraction`))
  }
}

projected <- PACEProjector(beta, proportionOfProbesRequired=opts$`min-probe-fraction`)
if (!is.list(projected) || length(projected) == 0) {
  stop('PACEProjector did not return its expected named model list')
}

model_names <- names(projected)
preferred <- which(grepl('Age45|PACE', model_names, ignore.case=TRUE))
idx <- if (length(preferred) > 0) preferred[1] else 1
model_name <- model_names[idx]
if (is.null(model_name) || !nzchar(model_name)) stop('Selected DunedinPACE model has no name')
scores <- projected[[idx]]

if (is.null(mPACE_Models$model_intercept[[model_name]])) {
  stop(sprintf('No model intercept found for %s', model_name))
}
model_intercept <- as.numeric(mPACE_Models$model_intercept[[model_name]])
if (length(model_intercept) != 1 || !is.finite(model_intercept)) {
  stop(sprintf('Invalid model intercept for %s', model_name))
}

parse_meta <- function(id) {
  m <- regexec('^(WT|SIRT[1-7])_rep([12])', id)
  z <- regmatches(id, m)[[1]]
  if (length(z) != 3) stop(paste('sample id must start CONDITION_repN:', id))
  list(condition=z[2], replicate=as.integer(z[3]))
}

sample_names <- names(scores)
if (is.null(sample_names)) sample_names <- colnames(beta)
meta <- lapply(sample_names, parse_meta)
out <- data.table(
  sample_id=sample_names,
  condition=vapply(meta, `[[`, '', 'condition'),
  replicate=vapply(meta, `[[`, integer(1), 'replicate'),
  dunedinpace=as.numeric(scores),
  model=model_name
)
if (any(!is.finite(out$dunedinpace))) stop('PACEProjector returned non-finite scores')

model_meta <- data.table(
  model=model_name,
  intercept=model_intercept,
  required_background_probes=length(required),
  annotation='IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
  source_package='danbelsky/DunedinPACE'
)

fwrite(out, opts$scores)
fwrite(rbindlist(qc), opts$qc)
fwrite(model_meta, opts$`model-metadata`)
