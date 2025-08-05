#!/usr/bin/env Rscript

# scripts/seurat_to_anndata_F.R
# Convert female Seurat objects to AnnData format for Python analysis

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

`%!in%` <- Negate(`%in%`)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample ID as argument")
}

sample_id <- args[1]
filt_path <- 'mapped_reads_F'
out_path <- 'h5ad'

cat("Converting female sample:", sample_id, "to AnnData format\n")

# Create output directory
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# Load filtered Seurat object
input_file <- file.path(filt_path, sample_id, paste0(sample_id, "_bender_filtered.rds"))
if (!file.exists(input_file)) {
  stop("Filtered Seurat object not found: ", input_file)
}

cat("Loading filtered Seurat object...\n")
loaded <- readRDS(input_file)

# Extract metadata and counts
cat("Extracting metadata and counts...\n")
m <- loaded@meta.data
c <- loaded[["RNA"]]$counts

# Create new Seurat object with clean structure
cat("Creating clean Seurat object...\n")
loadednow <- CreateSeuratObject(counts = c, meta.data = m)
DefaultAssay(loadednow) <- 'RNA'
loadednow[["RNA"]] <- as(object = loadednow[["RNA"]], Class = "Assay")

# Save as h5Seurat
h5seurat_file <- file.path(out_path, paste0(sample_id, "_bender_filtered.h5Seurat"))
cat("Saving H5Seurat file:", h5seurat_file, "\n")
SaveH5Seurat(loadednow, filename = h5seurat_file, overwrite = TRUE)

# Convert to h5ad
h5ad_file <- file.path(out_path, paste0(sample_id, "_bender_filtered.h5ad"))
cat("Converting to AnnData format:", h5ad_file, "\n")
Convert(h5seurat_file, dest = "h5ad", overwrite = TRUE)

# Clean up intermediate file
if (file.exists(h5seurat_file)) {
  file.remove(h5seurat_file)
}

cat("Successfully converted sample:", sample_id, "to AnnData format\n")
cat("Output file:", h5ad_file, "\n")
