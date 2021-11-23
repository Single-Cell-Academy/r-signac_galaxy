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
    help = "Filtered peak BC matrix file in h5 format."
  ),
  make_option(
    c("--ens-db-genome"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Metadata file."
  ),
  make_option(
    c("--annotations"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Genome version."
  ),
  make_option(
    c("--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R matrix object."
  ),
  make_option(
    c("--annotations-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Annotations file"
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
sessionInfo()
signac_object <- readRDS(file = opt$signac_object)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "GRCh38"

# add the gene information to the object
Annotation(signac_object) <- annotations

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
