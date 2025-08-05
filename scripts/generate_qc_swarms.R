#!/usr/bin/env Rscript

# scripts/generate_qc_swarms.R
# Generate swarm files for QC and doublet removal

library(stringr)

cat("Generating swarm files for QC and doublet removal\n")

# Generate swarm files for male samples
if (file.exists('sample_list_M.txt')) {
  all_files <- read.table('sample_list_M.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/qc-M.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "qc-M.swarm")
  cat("Created qc-M.swarm with", length(swarm_lines), "jobs\n")
}

# Generate swarm files for female samples
if (file.exists('sample_list_F.txt')) {
  all_files <- read.table('sample_list_F.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/qc-F.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "qc-F.swarm")
  cat("Created qc-F.swarm with", length(swarm_lines), "jobs\n")
}

cat("QC swarm file generation completed\n")
