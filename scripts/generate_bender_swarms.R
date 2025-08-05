#!/usr/bin/env Rscript

# scripts/generate_bender_swarms.R
# Generate swarm files for CellBender to Seurat conversion

library(stringr)

cat("Generating swarm files for CellBender to Seurat conversion\n")

# Generate swarm files for male samples
if (file.exists('sample_list_M.txt')) {
  all_files <- read.table('sample_list_M.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/bender-seu-M.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "bender-seu-M.swarm")
  cat("Created bender-seu-M.swarm with", length(swarm_lines), "jobs\n")
}

# Generate swarm files for female samples
if (file.exists('sample_list_F.txt')) {
  all_files <- read.table('sample_list_F.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/bender-seu-F.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "bender-seu-F.swarm")
  cat("Created bender-seu-F.swarm with", length(swarm_lines), "jobs\n")
}

cat("CellBender to Seurat swarm file generation completed\n")
