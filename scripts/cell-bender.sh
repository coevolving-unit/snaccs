#!/bin/bash

# scripts/cell-bender.sh
# Convert SCE objects to h5ad format and run CellBender for ambient RNA removal

set -e

echo "Starting CellBender ambient RNA removal pipeline"

# Step 1: Convert SingleCellExperiment objects to h5ad format
echo "Step 1: Converting SCE objects to h5ad format"

# Create conversion scripts
cat > scripts/sce-to-h5ad-M.R << 'EOF'
#!/usr/bin/env Rscript

# scripts/sce-to-h5ad-M.R
# Convert male SCE objects to h5ad format for CellBender

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
filt_path <- 'mapped_reads_M'

cat("Converting male sample:", sample_id, "to h5ad format\n")

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
EOF

cat > scripts/sce-to-h5ad-F.R << 'EOF'
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
EOF

# Create swarm file generation script
cat > scripts/generate_conversion_swarms.R << 'EOF'
#!/usr/bin/env Rscript

# scripts/generate_conversion_swarms.R
# Generate swarm files for SCE to h5ad conversion

library(stringr)

cat("Generating swarm files for SCE to h5ad conversion\n")

# Generate swarm files for male samples
if (file.exists('sample_list_M.txt')) {
  all_files <- read.table('sample_list_M.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/sce-to-h5ad-M.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "sce-to-h5ad-M.swarm")
  cat("Created sce-to-h5ad-M.swarm with", length(swarm_lines), "jobs\n")
}

# Generate swarm files for female samples
if (file.exists('sample_list_F.txt')) {
  all_files <- read.table('sample_list_F.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/sce-to-h5ad-F.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "sce-to-h5ad-F.swarm")
  cat("Created sce-to-h5ad-F.swarm with", length(swarm_lines), "jobs\n")
}

cat("Swarm file generation completed\n")
EOF

# Generate swarm files
Rscript scripts/generate_conversion_swarms.R

# Submit conversion jobs
echo "Submitting SCE to h5ad conversion jobs"
swarm -f sce-to-h5ad-M.swarm --module R/4.3 -g 10
swarm -f sce-to-h5ad-F.swarm --module R/4.3 -g 10

# Wait for conversion jobs to complete
echo "Waiting for conversion jobs to complete..."
echo "Please check job status with 'squeue -u \$USER' before proceeding to Step 2"
echo ""
echo "Once conversion is complete, run Step 2:"
echo "bash scripts/run_cellbender_step2.sh"
