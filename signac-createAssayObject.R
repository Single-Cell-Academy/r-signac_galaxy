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
    help = "A Seurat object."
  ),
  make_option(
    c("--counts"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Counts Matrix."
  ),
  make_option(
    c("--name"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Assay Name."
  ),
  make_option(
    c("--min-cells"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Min Cells."
  ),
  make_option(
    c("--min-features"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Min Features."
  ),
  make_option(
    c("--method"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Method for normalization."
  ),
  make_option(
    c("--scale-factor"),
    action = "store",
    default = 10000,
    type = 'numeric',
    help = "Sets the scale factor for cell-level normalization."
  ),
  make_option(
    c("--margin"),
    action = "store",
    default = NA,
    type = 'character',
    help = "If performing CLR normalization, normalize across features (1) or cells (2)."
  ),
  make_option(
    c("--output-object-file"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "File name in which to store serialized R matrix object."
  ) 
)
opt <- wsc_parse_args(option_list)

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(Signac))

set.seed(1234)

if (! file.exists(opt$signac_object)){
  stop((paste('File', opt$signac_object, 'does not exist')))
}
if (! file.exists(opt$counts)){
  stop((paste('File', opt$counts, 'does not exist')))
}

signac_object <- readRDS(file = opt$signac_object)
counts <- readRDS(opt$counts)

signac_object[[opt$name]] <- CreateAssayObject(counts = counts)#, min.cells = opt$min_cells, min.features = opt$min_features)

signac_object <- NormalizeData(object = signac_object, normalization.method = opt$method, scale.factor = opt$scale_factor, margin = opt$margin, assay = opt$name)

saveRDS(signac_object, file = opt$output_object_file)
