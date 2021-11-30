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
    c("--reduction-use"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--dims-use"),
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

dims_use <- opt$dims_use
if ( ! is.null(dims_use)){
  dims_use <- as.integer(wsc_parse_numeric(opt, 'dims_use'))
}

suppressPackageStartupMessages(require(Signac))

set.seed(1234)

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)

signac_object <- RunUMAP(object = signac_object, reduction = opt$reduction_use, dims = dims_use)

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
