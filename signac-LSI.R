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
    c("--min-cutoff"),
    action = "store",
    default = NA,
    type = 'character',
    help = ""
  ),
  make_option(
    c("--output-depthcor"),
    action = "store",
    default = NA,
    type = 'character',
    help = "TSS output plot."
  ),
  make_option(
    c("-w", "--png-width"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Width of png (px)."
  ),
  make_option(
    c("-j", "--png-height"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Height of png file (px)."
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

set.seed(1234)

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)


signac_object <- RunTFIDF(signac_object)
signac_object <- FindTopFeatures(signac_object, min.cutoff = opt$min_cutoff)
signac_object <- RunSVD(signac_object)

## Plot the Depth correlation plot
png(filename = opt$output_depthcor, width = opt$png_width, height = opt$png_height)
DepthCor(signac_object)
dev.off()

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
