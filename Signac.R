rm(list = ls())

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

library(future)
plan("multicore", workers = (availableCores()-1))
options(future.globals.maxSize = 3000000 * 1024^2)

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/scATACseq/10XData/atac_v1_pbmc_10k_output/")

counts <- Read10X_h5(filename = "./atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")


metadata <- read.csv(
  file = "./atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

metadata <- metadata[colnames(counts),]

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'GRCh38',
  fragments = './atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/scATACseq/Signac_analysis/atac_pbmc_500_nextgem")

pbmc[['peaks']]

granges(pbmc)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to GRCh38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "GRCh38"

# add the gene information to the object
Annotation(pbmc) <- annotations


# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments



pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')

png("Tssplot.png")
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
dev.off()



pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
png("FragmentHistogram.png")
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
dev.off()

png("VlnPlot_QC.png", width=1000)
VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()


pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

png("DepthCor.png")
DepthCor(pbmc)
dev.off()

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)


png("UMAP.png")
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()


gene.activities <- GeneActivity(pbmc)


# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)


DefaultAssay(pbmc) <- 'RNA'

png("FeaturePlot_knownMarkers.png", width=1000)
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()


saveRDS(pbmc, paste("Seurat_object.rds", sep="."))


















q("no")
