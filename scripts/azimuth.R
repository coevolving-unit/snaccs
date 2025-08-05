#!/usr/bin/env Rscript

# azimuth.R
# Annotate input Seurat object using a reference with Azimuth

# Load necessary libraries
library(Azimuth)
library(Seurat)
library(anndata)
library(SeuratDisk)
library(Matrix)
library(SeuratData)
library(ggplot2)
library(stringr)
library(BPCells)
library(scrattch.io)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
class <- args[1]
ref <- args[2]

cat("Loading data...\n")

# Read h5ad and metadata
f <- read_h5ad(paste0(class, "_annotate.h5ad"), backed = NULL)
f$X <- Matrix(f$X, sparse = TRUE)

m <- read.csv(paste0(class, "_annotate_meta.csv"))
rownames(m) <- m$X

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = t(f$X), meta.data = m)
str(seurat_obj)

cat("Running Azimuth...\n")

# Run Azimuth
seurat_obj_annotate <- RunAzimuth(seurat_obj, reference = ref)

cat("Saving output...\n")

# Save annotated metadata
saveRDS(seurat_obj_annotate@meta.data, file = paste0(class, "_", ref, ".rds"))

cat("Done.\n")
