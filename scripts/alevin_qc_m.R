#!/usr/bin/env Rscript

# scripts/alevin_qc_m.R
# Generate QC reports for male samples

library('alevinQC')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample ID as argument")
}

sample_id <- args[1]
sexnow <- 'mapped_reads_M'

cat("Generating QC report for male sample:", sample_id, "\n")

# Define paths
map_dir <- paste(sexnow, paste0(sample_id, "_map"), sep = "/")
permit_dir <- paste(sexnow, paste0(sample_id, "_quant"), sep = "/")
quant_dir <- paste(sexnow, paste0(sample_id, "_quant_res"), sep = "/")

# Check if files exist
check <- checkAlevinFryInputFiles(map_dir, permit_dir, quant_dir)
print(check)

# Generate QC report
alevinFryQCReport(
  mapDir = map_dir,
  permitDir = permit_dir,
  quantDir = quant_dir,
  sampleId = sample_id,
  outputFormat = 'html_document',
  outputFile = paste0(sample_id, ".html"),
  outputDir = paste0(sexnow, "/"),
  forceOverwrite = TRUE
)

cat("QC report completed for sample:", sample_id, "\n")
