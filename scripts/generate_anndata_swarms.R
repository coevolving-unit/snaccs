#!/usr/bin/env Rscript

# scripts/generate_anndata_swarms.R
# Generate swarm files for Seurat to AnnData conversion

library(stringr)

cat("Generating swarm files for Seurat to AnnData conversion\n")

# Generate swarm files for male samples
if (file.exists('sample_list_M.txt')) {
  all_files <- read.table('sample_list_M.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/seurat_to_anndata_M.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "seurat_to_anndata_M.swarm")
  cat("Created seurat_to_anndata_M.swarm with", length(swarm_lines), "jobs\n")
}

# Generate swarm files for female samples
if (file.exists('sample_list_F.txt')) {
  all_files <- read.table('sample_list_F.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/seurat_to_anndata_F.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "seurat_to_anndata_F.swarm")
  cat("Created seurat_to_anndata_F.swarm with", length(swarm_lines), "jobs\n")
}

cat("AnnData conversion swarm file generation completed\n")
