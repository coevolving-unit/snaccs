#!/usr/bin/env Rscript

# scripts/sce-to-h5ad-F.R  
# Convert female SCE objects to h5ad format for CellBender

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample ID as argument")
}

sample_id <- args[1]
filt_path <- 'mapped_reads_F'

cat("Converting female sample:", sample_id, "to h5ad format\n")

# Load SCE object
sce_file <- file.path(filt_path, paste0(sample_id, "_sce.rds"))
if (!file.exists(sce_file)) {
  stop("SCE file not found: ", sce_file)
}

loaded <- readRDS(sce_file)

# Convert to Seurat object
loaded <- CreateSeuratObject(counts = round(counts(loaded), digits = 0))
loaded[["RNA"]] <- as(object = loaded[["RNA"]], Class = "Assay")

# Save as h5Seurat
h5seurat_file <- file.path(filt_path, paste0(sample_id, ".h5Seurat"))
SaveH5Seurat(loaded, filename = h5seurat_file, overwrite = TRUE)

# Convert to h5ad
Convert(h5seurat_file, dest = "h5ad", overwrite = TRUE)

cat("Conversion completed for sample:", sample_id, "\n")
