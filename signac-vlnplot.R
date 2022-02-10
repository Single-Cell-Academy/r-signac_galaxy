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
    c("--features"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File to be used to derive a vector of genes to test."
  ),
  make_option(
    c("--cols"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of colors."
  ),
  make_option(
    c("-p", "--pt-size"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Adjust point size for plotting."
  ),
  make_option(
    c("--group-by"),
    action = "store",
    default = "ident",
    type = 'character',
    help = "Group (color) cells in different ways (for example, orig.ident)."
  ),
  make_option(
    c("--sort"),
    action = "store",
    type = 'logical',
    help = "Sort identity classes (on the x-axis) by the average expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction."
  ),
  make_option(
    c("--assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Name of assay to use, defaults to the active assay."
  ),
  # make_option(
  #   c("--split-by"),
  #   action = "store",
  #   default = NULL,
  #   type = 'character',
  #   help = "A variable to split the violin plots by."
  # ),
  make_option(
    c("--same-y-lims"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Set all the y-axis limits to the same values."
  ),
  make_option(
    c("--log"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "plot the feature axis on log scale."
  ),
  make_option(
    c("--ncol"),
    action = "store",
    default = NULL,
    type = 'integer',
    help = "Number of columns if multiple plots are displayed."
  ),
  make_option(
    c("--slot"),
    action = "store",
    default = 'data',
    type = 'character',
    help = "Use non-normalized counts data for plotting."
  ),
  # make_option(
  #   c("--split-plot"),
  #   action = "store",
  #   default = FALSE,
  #   type = 'logical',
  #   help = "plot each group of the split violin plots by multiple or single violin shapes."
  # ),
  make_option(
    c("--stack"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Horizontally stack plots for each feature."
  ),
  make_option(
    c("--fill-by"),
    action = "store",
    default = 'feature',
    type = 'character',
    help = "Color violins/ridges based on either 'feature' or 'ident'."
  ),  
  make_option(
    c("--flip"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "flip plot orientation (identities on x-axis)."
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

if (! file.exists(opt$signac_object)){
  stop((paste('File', opt$signac_object, 'does not exist')))
}

# Check features
features <- NULL
if (! is.null(opt$features) && opt$features != 'NULL'){
  if (! file.exists(opt$features)){
    stop((paste('Supplied features file', opt$features, 'does not exist')))
  }else{
    features <- readLines(opt$features)
  }
}

# Check cols
cols <- NULL
if (! is.null(opt$cols) && opt$cols != 'NULL'){
  if (! file.exists(opt$cols)){
    stop((paste('Supplied features file', opt$cols, 'does not exist')))
  }else{
    cols <- readLines(opt$cols)
  }
}

# Check assay
assay <- NULL
if (! is.null(opt$assay) && opt$assay != 'NULL'){
  assay <- opt$assay
}

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)

group_by <- if(opt$group_by == 'ident'){
  NULL
}else{
  opt$group_by
}

# sort <- if(opt$sort == 'False'){
#   FALSE
# }else{
#   TRUE
# }

png(filename = opt$output_image_file, width = opt$png_width, height = opt$png_height)
VlnPlot(object = signac_object, 
  features = features, 
  cols = opt$cols, 
  pt.size = opt$pt_size, 
  group.by = group_by,
  sort = opt$sort,
  assay = assay, 
  #split.by = opt$split_by,
  same.y.lims = opt$same_y_lims, 
  log = opt$log, 
  ncol = opt$ncol,
  slot = opt$slot, 
  #split.plot = opt$split_plot, 
  stack = opt$stack,
  fill.by = opt$fill_by, 
  flip = opt$flip
  )
dev.off()