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
    c("--"),
    action = "store",
    default = NA,
    type = 'character',
    help = "."
  ),
  make_option(
    c("--"),
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



# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
