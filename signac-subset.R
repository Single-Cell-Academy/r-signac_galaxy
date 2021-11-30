#!/usr/bin/env Rscript

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("--signac-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = ""
  ),
  make_option(
    c("--peak-region-fragments-min"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--peak-region-fragments-max"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--pct-reads-in-peaks"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--blacklist-ratio"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--nucleosome-signal"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--TSS-enrichment"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R matrix object."
  ) 
)

opt <- wsc_parse_args(option_list)

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(Signac))
suppressPackageStartupMessages(require(GenomeInfoDb))
suppressPackageStartupMessages(require(EnsDb.Hsapiens.v75))
suppressPackageStartupMessages(require(EnsDb.Mmusculus.v79))

set.seed(1234)

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)

signac_object <- subset(
  x = signac_object,
  subset = peak_region_fragments > opt$peak_region_fragments_min &
    peak_region_fragments < opt$peak_region_fragments_max &
    pct_reads_in_peaks > opt$pct_reads_in_peaks &
    blacklist_ratio < opt$blacklist_ratio &
    nucleosome_signal < opt$nucleosome_signal &
    TSS.enrichment > opt$TSS_enrichment
)

signac_object

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
