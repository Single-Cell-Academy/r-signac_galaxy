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
    c("-a", "--dim-1"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Dimension for x-axis (default 1)"
  ),
make_option(
    c("-b", "--dim-2"),
    action = "store",
    default = 2,
    type = 'integer',
    help = "Dimension for y-axis (default 2)"
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
    c("-f", "--group-by"),
    action = "store",
    default = 'ident',
    type = 'character',
    help = "Group (color) cells in different ways (for example, orig.ident)."
  ),
make_option(
    c("-t", "--plot-title"),
    action = "store",
    default = 1,
    type = 'character',
    help = "Title for plot."
  ),
make_option(
    c("-m", "--do-bare"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Do only minimal formatting (default : FALSE)"
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

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)

png(filename = opt$output_image_file, width = opt$png_width, height = opt$png_height)
DimPlot(object = signac_object, dim.1 = opt$dim_1, dim.2 = opt$dim_2, pt.size = opt$pt_size, label.size = opt$label_size, do.label = opt$do_label, group.by = opt$group_by, do.bare=opt$do_bare) + NoLegend()
dev.off()
