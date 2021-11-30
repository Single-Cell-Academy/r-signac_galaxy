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
    c("--tss-enrichment"),
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

suppressPackageStartupMessages(require(Signac))

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)


## transform input parameters to numeric
peak_region_fragments_min <- as.numeric(opt$peak_region_fragments_min)
peak_region_fragments_max <- as.numeric(opt$peak_region_fragments_max)
pct_reads_in_peaks_var <- as.numeric(opt$pct_reads_in_peaks)
blacklist_ratio_var <- as.numeric(opt$blacklist_ratio)
nucleosome_signal_var <- as.numeric(opt$nucleosome_signal)
tss_enrichment_var <- as.numeric(opt$tss_enrichment)

print("Signac object before filtering:")

signac_object

signac_object <- subset(signac_object, peak_region_fragments > peak_region_fragments_min)
signac_object <- subset(signac_object, peak_region_fragments < peak_region_fragments_max)
signac_object <- subset(signac_object, pct_reads_in_peaks > pct_reads_in_peaks_var)
signac_object <- subset(signac_object, blacklist_ratio < blacklist_ratio_var)
signac_object <- subset(signac_object, nucleosome_signal < nucleosome_signal_var)
signac_object <- subset(signac_object, TSS.enrichment > tss_enrichment_var)

print("Signac object after filtering:")

signac_object

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
