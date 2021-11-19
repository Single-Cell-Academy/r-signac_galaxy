#!/usr/bin/env Rscript

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("--h5-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Filtered peak BC matrix file in h5 format."
  ),
  make_option(
    c("--metadata"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Metadata file."
  ),
  make_option(
    c("--genome"),
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
    c("--fragment-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Fragment file from CellRanger-ATAC."
  )
)

opt <- wsc_parse_args(option_list)

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(Signac))

atac_h5_matrix <- Read10X_h5(filename = opt$h5_file)

metadata <- read.csv(
  file = opt$metadata,
  header = TRUE,
  row.names = 1
)
metadata <- subset(metadata,rownames(metadata) %in% colnames(atac_h5_matrix))

chrom_assay <- CreateChromatinAssay(
  counts = atac_h5_matrix,
  sep = c(":", "-"),
  genome = opt$genome,
  fragments = opt$fragment_file,
  min.cells = 0,
  min.features = 0
)

signac_object <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# cat(c(
#   '# Object summary', 
#   capture.output(print(signac_object)), 
#   '\n# Metadata sample', 
#   capture.output(head(signac_object@meta.data))
# ), 
# sep = '\n')

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
