#!/usr/bin/env Rscript

# scripts/bender-seu-M.R
# Load CellBender h5 files and convert to Seurat objects for males
# Includes Ensembl ID to gene name conversion and % mito calculation

suppressPackageStartupMessages({
  library(scCustomize)
  library(Seurat)
  library(SeuratDisk)
  library(SingleCellExperiment)
  library(dplyr)
})

`%!in%` <- Negate(`%in%`)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample ID as argument")
}

sample_id <- args[1]
filt_path <- 'mapped_reads_M'

cat("Processing male sample:", sample_id, "\n")

# Load Cell Ranger gene reference for ID conversion
# Use output from cellranger.sh script
cr_output_dir <- "cellranger_output"
cr_genes_file <- NULL

# Try to find any Cell Ranger output to use as gene reference
if (dir.exists(cr_output_dir)) {
  # Look for any sample's features file to use as gene reference
  sample_dirs <- list.dirs(cr_output_dir, recursive = FALSE)
  for (sample_dir in sample_dirs) {
    potential_file <- file.path(sample_dir, "outs/filtered_feature_bc_matrix/features.tsv.gz")
    if (file.exists(potential_file)) {
      cr_genes_file <- potential_file
      break
    }
  }
}

# Fallback to a specific sample if available
if (is.null(cr_genes_file)) {
  # Try using the current sample's Cell Ranger output
  sample_cr_file <- file.path(cr_output_dir, sample_id, "outs/filtered_feature_bc_matrix/features.tsv.gz")
  if (file.exists(sample_cr_file)) {
    cr_genes_file <- sample_cr_file
  }
}

# Final fallback to original hardcoded path if needed
if (is.null(cr_genes_file)) {
  cr_path <- "/data/NIMH_scratch/GTEx/process/snRNAseq_cellranger/Jenny04282022/TRIMMED"
  cr_genes_file <- file.path(cr_path, "sample_80914/outs/filtered_feature_bc_matrix/features.tsv.gz")
  cat("Warning: Using fallback gene reference path\n")
}

if (!file.exists(cr_genes_file)) {
  stop("Cell Ranger genes file not found: ", cr_genes_file, 
       "\nPlease ensure Cell Ranger has been run or features.tsv.gz is available")
}

crgenes <- read.table(cr_genes_file, stringsAsFactors = FALSE)

cat("Loading and converting original dataset\n")

# Load original SCE object
ori_file <- file.path(filt_path, paste0(sample_id, "_sce.rds"))
if (!file.exists(ori_file)) {
  stop("Original SCE file not found: ", ori_file)
}

ori <- readRDS(ori_file)
gn <- as.matrix(counts(ori))

# Convert Ensembl IDs to gene names
idx <- match(rownames(gn), crgenes$V1)
gi <- crgenes[idx, ]

# Handle duplicate gene names
nn <- gi %>% group_by(V2) %>% summarise(n = n())
nn <- subset(nn, n == 2)
gs <- subset(gi, V2 %in% nn$V2)

gn1 <- gn[which(rownames(gn) %in% gs$V1), ]
gn2 <- gn[which(rownames(gn) %!in% gs$V1), ]
na <- subset(gi, V1 %!in% gs$V1)

gs <- gs[match(rownames(gn1), gs$V1), ]
na <- na[match(rownames(gn2), na$V1), ]
rownames(gn2) <- na$V2

rm(gn)

# Aggregate duplicate gene names
gn1b <- Matrix.utils::aggregate.Matrix(gn1, groupings = gs$V2, fun = "sum")
gn <- rbind(gn1b, gn2)

# Create Seurat object
ori <- CreateSeuratObject(counts = gn)
rm(gn)

# Calculate mitochondrial percentage
ori[["percent.mt"]] <- PercentageFeatureSet(object = ori, pattern = "^MT-")

ori <- ori@meta.data
cat("Original dataset dimensions:", nrow(ori), "x", ncol(ori), "\n")
cat("Cells with >5% mito:", sum(ori$percent.mt > 5), "/", nrow(ori), "\n")

cat("Loading and converting filtered CellBender dataset\n")

# Load CellBender filtered results
bender_file <- file.path(filt_path, sample_id, paste0(sample_id, "_bender_filtered.h5"))
if (!file.exists(bender_file)) {
  # Try alternative naming from cellbender_cuda output
  bender_file <- file.path(filt_path, sample_id, paste0(sample_id, "_bender.h5"))
  if (!file.exists(bender_file)) {
    stop("CellBender file not found. Tried: ", 
         file.path(filt_path, sample_id, paste0(sample_id, "_bender_filtered.h5")), " and ", 
         bender_file)
  }
}

loaded <- Read_CellBender_h5_Mat(file_name = bender_file)
gn <- as.matrix(loaded)

# Convert Ensembl IDs to gene names (same process as above)
idx <- match(rownames(gn), crgenes$V1)
gi <- crgenes[idx, ]
nn <- gi %>% group_by(V2) %>% summarise(n = n())
nn <- subset(nn, n == 2)
gs <- subset(gi, V2 %in% nn$V2)

gn1 <- gn[which(rownames(gn) %in% gs$V1), ]
gn2 <- gn[which(rownames(gn) %!in% gs$V1), ]
na <- subset(gi, V1 %!in% gs$V1)

if (nrow(gn2) > 0 && nrow(na) > 0) {
  rownames(gn2) <- na$V2
}

rm(gn)

gn1b <- Matrix.utils::aggregate.Matrix(gn1, groupings = gs$V2, fun = "sum")
gn1b <- as.matrix(gn1b)
gnn <- rbind(gn1b, gn2)

# Create filtered Seurat object
loaded <- CreateSeuratObject(counts = gnn)
loaded[["percent.mt"]] <- PercentageFeatureSet(object = loaded, pattern = "^MT-")

cat("Filtered dataset dimensions:", nrow(loaded@meta.data), "x", ncol(loaded@meta.data), "\n")
cat("Cells with >5% mito:", sum(loaded@meta.data$percent.mt > 5), "/", nrow(loaded@meta.data), "\n")

cat("Adding original counts/features to CellBender dataset\n")

# Add original metrics to filtered dataset
ori <- ori[which(rownames(ori) %in% rownames(loaded@meta.data)), ]
idx <- match(rownames(loaded@meta.data), rownames(ori))
ori <- ori[idx, ]

loaded@meta.data$nCount_RNA_ori <- ori$nCount_RNA
loaded@meta.data$nFeature_RNA_ori <- ori$nFeature_RNA
loaded@meta.data$percent.mt_ori <- ori$percent.mt

# Save processed object
output_file <- file.path(filt_path, sample_id, paste0(sample_id, "_bender.rds"))
saveRDS(loaded, output_file)

cat("Saved processed object to:", output_file, "\n")
cat("Processing completed for sample:", sample_id, "\n")
