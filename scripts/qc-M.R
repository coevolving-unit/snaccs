#!/usr/bin/env Rscript

# scripts/qc-M.R
# Remove cells > 5% mito, remove doublets, and make QC plots for males

suppressPackageStartupMessages({
  library(Seurat)
  library(scater)
  library(SingleCellExperiment)
  library(scDblFinder)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample ID as argument")
}

sample_id <- args[1]
filt_path <- 'mapped_reads_M'

cat("Processing QC for male sample:", sample_id, "\n")

# Load processed Seurat object
input_file <- file.path(filt_path, sample_id, paste0(sample_id, "_bender.rds"))
if (!file.exists(input_file)) {
  stop("Bender processed file not found: ", input_file)
}

loaded <- readRDS(input_file)

# Remove high mito cells
cat("Removing high mito cells (>5%)\n")
initial_cells <- ncol(loaded)
loaded <- subset(loaded, percent.mt < 5)
remaining_cells <- ncol(loaded)
cat("Removed", initial_cells - remaining_cells, "cells with >5% mito\n")
cat("Remaining cells:", remaining_cells, "\n")

# Normalize and process for clustering
cat("Converting to SCE and processing\n")

loaded <- NormalizeData(loaded)
loaded <- FindVariableFeatures(loaded, selection.method = "vst", nfeatures = 2000)
loaded <- ScaleData(loaded)
loaded <- RunPCA(loaded)

# Determine optimal number of PCs
pct <- loaded[["pca"]]@stdev / sum(loaded[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs <- min(co1, co2)

cat("Using", pcs, "principal components\n")

# Run dimensionality reduction and clustering
loaded <- RunTSNE(loaded, dims = 1:pcs, check_duplicates = FALSE)
loaded <- RunUMAP(loaded, dims = 1:pcs)
loaded <- FindNeighbors(loaded, dims = 1:pcs) 
loaded <- FindClusters(loaded, resolution = 1.2)

# Convert to SCE for doublet detection
sce <- as.SingleCellExperiment(loaded)
sce <- logNormCounts(sce)

cat("Running doublet detection\n")

# Calculate expected doublet rate
exp_rate <- (0.0527 + 0.0008 * ncol(loaded)) / 100
cat("Expected doublet rate:", round(exp_rate * 100, 2), "%\n")

# Run scDblFinder
sce.dbl <- scDblFinder(sce, clusters = sce@colData$ident, dbr = exp_rate)
loaded@meta.data <- data.frame(sce.dbl@colData)

# Count doublets
doublet_count <- sum(loaded@meta.data$scDblFinder.class == "doublet")
cat("Detected doublets:", doublet_count, "/", ncol(loaded), "\n")

cat("Creating QC plots\n")

# Create QC plots
pdf_file <- file.path(filt_path, paste0(sample_id, "_QC.pdf"))
pdf(file = pdf_file, width = 12, height = 8)

# Clustering and doublet plots
DimPlot(loaded, reduction = 'umap', group.by = "RNA_snn_res.1.2")
DimPlot(loaded, reduction = 'umap', group.by = "scDblFinder.class")
FeaturePlot(loaded, reduction = 'umap', features = c("scDblFinder.score"), label = TRUE)

# QC metric plots
VlnPlot(loaded, features = "percent.mt", split.by = "RNA_snn_res.1.2")
FeaturePlot(loaded, reduction = 'umap', features = c("percent.mt"), label = TRUE)
VlnPlot(loaded, features = "nCount_RNA", split.by = "RNA_snn_res.1.2")
FeaturePlot(loaded, reduction = 'umap', features = c("nCount_RNA"), label = TRUE)
VlnPlot(loaded, features = "intronRat", split.by = "RNA_snn_res.1.2")
FeaturePlot(loaded, reduction = 'umap', features = c("intronRat"), label = TRUE)

# Ambient RNA markers
ambient_genes <- c("MALAT1", "KCNIP4", "CSMD1", "RBFOX1", "RALYL", "SYT1", "NRGN", "CHN1")
for (gene in ambient_genes) {
  if (gene %in% rownames(loaded)) {
    FeaturePlot(loaded, reduction = 'umap', features = gene, label = TRUE)
  }
}

# Cell type markers
marker_genes <- list(
  oligo = c("OLIG1", "OLIG2", "PLP1"),
  opc = c("SOX10", "CSPG4"),
  astro = c("GFAP", "ALDH1L1", "AQP4"),
  micro = c("CD74", "APBB1IP"),
  immune = c("MRC1", "PTPRC"),
  neurons = c("NEUROD6", "GRIN1", "GRIN2B", "SLC17A7", "SLC17A6"),
  gaba = c("GAD1", "CALB2", "LAMP5", "SST", "VIP", "PVALB"),
  vascular = c("MCAM", "FLT1", "DCN")
)

for (cell_type in names(marker_genes)) {
  for (gene in marker_genes[[cell_type]]) {
    if (gene %in% rownames(loaded)) {
      FeaturePlot(loaded, reduction = 'umap', features = gene, label = TRUE)
    }
  }
}

dev.off()

cat("QC plots saved to:", pdf_file, "\n")

cat("Removing doublet cells\n")
loaded <- subset(loaded, scDblFinder.class == 'singlet')
final_cells <- ncol(loaded)
cat("Final cell count after doublet removal:", final_cells, "\n")

# Save filtered object
output_file <- file.path(filt_path, sample_id, paste0(sample_id, "_bender_filtered.rds"))
saveRDS(loaded, output_file)

cat("Saved filtered object to:", output_file, "\n")
cat("QC processing completed for sample:", sample_id, "\n")
