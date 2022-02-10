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
    c("--fragment-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Fragments file."
  ),
  make_option(
    c("--assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Assay."
  ),
  make_option(
    c("--features"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Features."
  ),
  make_option(
    c("--extend-upstream"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Number of bases to extend upstream of the TSS."
  ),
  make_option(
    c("--extend-downstream"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Number of bases to extend downstream of the TSS."
  ),
  make_option(
    c("--biotypes"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Gene biotypes to include. If NULL, use all biotypes in the gene annotation."
  ),
  make_option(
    c("--max-width"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Maximum allowed gene width for a gene to be quantified."
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

set.seed(1234)

if (! file.exists(opt$signac_object)){
  stop((paste('File', opt$signac_object, 'does not exist')))
}

signac_object <- readRDS(file = opt$signac_object)

signac_object@assays$peaks@fragments[[1]]@path <- opt$fragment_file

# Check features
features <- NULL
if (! is.null(opt$features) && opt$features != 'NULL'){
  if (file.exists(opt$features)){
    features <- readLines(opt$features)
  }
}

# Check assay
assay <- NULL
if (! is.null(opt$assay) && opt$assay != 'NULL'){
  assay <- opt$assay
}

saveRDS(GeneActivity(object = signac_object, assay = assay, features = features, extend.upstream = opt$extend_upstream, extend.downstream = opt$extend_downstream, biotypes = opt$biotypes, max.width = opt$max_width), file = opt$output_object_file)
