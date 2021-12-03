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
    c("--dims"),
    action = "store",
    type = 'character',
    help = "Dimension for x-axis (default 1)"
  ),
  make_option(
    c("-p", "--pt-size"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Adjust point size for plotting"
  ),
  make_option(
    c("-l", "--label-size"),
    action = "store",
    default = 4,
    type = 'integer',
    help = "Sets size of labels."
  ),
  make_option(
    c("-d", "--do-label"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Whether to label the clusters."
  ),
  make_option(
    c("--group-by"),
    action = "store",
    default = 'ident',
    type = 'character',
    help = "Group (color) cells in different ways (for example, orig.ident)."
  ),
  make_option(
    c("--png-width"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Width of png (px)."
  ),
  make_option(
    c("--png-height"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Height of png file (px)."
  ),
  make_option(
    c("--output_image_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R matrix object."
  ) 
)

opt <- wsc_parse_args(option_list)

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(Signac))

set.seed(1234)

dims_use <- opt$dims
if ( ! is.null(dims_use)){
  dims_parsed <- wsc_parse_numeric(opt, 'dims')
  dims_use <- seq(from = dims_parsed[1], to = dims_parsed[2])
}

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)

png(filename = opt$output_image_file, width = opt$png_width, height = opt$png_height)
DimPlot(object = signac_object, dims = dims_use, pt.size = opt$pt_size, label.size = opt$label_size, group.by = opt$group_by) + NoLegend()
dev.off()
