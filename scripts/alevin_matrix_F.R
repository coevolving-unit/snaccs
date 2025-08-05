#!/usr/bin/env Rscript

# scripts/alevin_matrix_F.R
# Create SingleCellExperiment objects for female samples

# Source the load_fry function
source("scripts/utils/load_fry.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample ID as argument")
}

sample_id <- args[1]
cat("Processing sample:", sample_id, "\n")

# Set paths
quant_res <- paste0(sample_id, "_quant_res")
fry_path <- file.path("mapped_reads_F", quant_res)

# Load data
af <- load_fry(fry_path)

# Save SingleCellExperiment object
output_file <- file.path("mapped_reads_F", paste0(sample_id, "_sce.rds"))
saveRDS(af, file = output_file)

cat("Saved SCE object to:", output_file, "\n")
cat("Done processing sample:", sample_id, "\n")
