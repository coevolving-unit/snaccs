#!/usr/bin/env Rscript

# scripts/generate_matrix_swarms.R
# Generate swarm files for creating count matrices

library(stringr)

cat("Generating swarm files for matrix creation\n")

# Generate swarm files for male samples
if (file.exists('sample_list_M.txt')) {
  all_files <- read.table('sample_list_M.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/alevin_matrix_M.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "alevin_matrix_M_R.swarm")
  cat("Created alevin_matrix_M_R.swarm with", length(swarm_lines), "jobs\n")
}

# Generate swarm files for female samples
if (file.exists('sample_list_F.txt')) {
  all_files <- read.table('sample_list_F.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/alevin_matrix_F.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "alevin_matrix_F_R.swarm")
  cat("Created alevin_matrix_F_R.swarm with", length(swarm_lines), "jobs\n")
}

# Generate QC swarm files for male samples
if (file.exists('sample_list_M.txt')) {
  all_files <- read.table('sample_list_M.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/alevin_qc_m.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "alevin_qc_m.swarm")
  cat("Created alevin_qc_m.swarm with", length(swarm_lines), "jobs\n")
}

# Generate QC swarm files for female samples
if (file.exists('sample_list_F.txt')) {
  all_files <- read.table('sample_list_F.txt', stringsAsFactors = FALSE)
  swarm_lines <- character(0)
  
  for (i in 1:nrow(all_files)) {
    line <- paste("Rscript scripts/alevin_qc_f.R", all_files$V1[i])
    swarm_lines <- c(swarm_lines, line)
  }
  
  writeLines(swarm_lines, "alevin_qc_f.swarm")
  cat("Created alevin_qc_f.swarm with", length(swarm_lines), "jobs\n")
}

cat("Swarm file generation completed\n")
