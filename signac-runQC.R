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
    c("--tss_threshold"),
    action = "store",
    default = NA,
    type = 'character',
    help = "TSS enrichment threshold for marking regions as high tss regions."
  ),
  make_option(
    c("--output_tss_plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = "TSS output plot."
  ),
  make_option(
    c("--frag_history_plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Fragment length periodicity plot."
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

# compute nucleosome signal score per cell
signac_object <- NucleosomeSignal(object = signac_object)

# compute TSS enrichment score per cell
signac_object <- TSSEnrichment(object = signac_object, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
signac_object$pct_reads_in_peaks <- signac_object$peak_region_fragments / pbmc$passed_filters * 100
signac_object$blacklist_ratio <- signac_object$blacklist_region_fragments / pbmc$peak_region_fragments

signac_object$high.tss <- ifelse(signac_object$TSS.enrichment > tss_threshold, 'High', 'Low')

png(filename = opt$output_tss_plot)
TSSPlot(signac_object, group.by = 'high.tss') + NoLegend()
dev.off()

signac_object$nucleosome_group <- ifelse(signac_object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

png(filename = opt$frag_history_plot)
FragmentHistogram(object = signac_object, group.by = 'nucleosome_group')
dev.off()

# Output to a serialized R object
saveRDS(signac_object, file = opt$output_object_file)
