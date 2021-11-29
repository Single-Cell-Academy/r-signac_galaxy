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
    c("--tss-threshold"),
    action = "store",
    default = NA,
    type = 'character',
    help = "TSS enrichment threshold for marking regions as high tss regions."
  ),
  make_option(
    c("--output-tss-plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = "TSS output plot."
  ),
  make_option(
    c("--frag-history-plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Fragment length periodicity plot."
  ),
  make_option(
    c("--fragment-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Fragment file."
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

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(Signac))
suppressPackageStartupMessages(require(GenomeInfoDb))
suppressPackageStartupMessages(require(EnsDb.Hsapiens.v75))
suppressPackageStartupMessages(require(EnsDb.Mmusculus.v79))

set.seed(1234)

# extract gene annotations from EnsDb
signac_object <- readRDS(file = opt$signac_object)

## Modify fragments file location
current_wd <- getwd()
new_framgent_file_loc <- paste(current_wd,"fragments.tsv.gz",sep = "/")
signac_object@assays$peaks@fragments[[1]]@path <- new_framgent_file_loc

# compute nucleosome signal score per cell
signac_object <- NucleosomeSignal(object = signac_object)

# compute TSS enrichment score per cell
signac_object <- TSSEnrichment(object = signac_object, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
signac_object$pct_reads_in_peaks <- signac_object$peak_region_fragments / signac_object$passed_filters * 100
signac_object$blacklist_ratio <- signac_object$blacklist_region_fragments / signac_object$peak_region_fragments

signac_object$high.tss <- ifelse(signac_object$TSS.enrichment > opt$tss_threshold, 'High', 'Low')

png(filename = opt$output_tss_plot, width = opt$png_width, height = opt$png_height)
TSSPlot(signac_object, group.by = 'high.tss') + NoLegend()
dev.off()

signac_object$nucleosome_group <- ifelse(signac_object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

png(filename = opt$frag_history_plot, width = opt$png_width, height = opt$png_height)
FragmentHistogram(object = signac_object, group.by = 'nucleosome_group')
dev.off()

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
