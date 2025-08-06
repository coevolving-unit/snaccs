#!/usr/bin/env Rscript

# Process cellsnp-lite results into summary format
# Usage: Rscript 06_process_cellsnp_results.R <SAMPLE_NAME>

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(data.table)
library(Matrix)

# Parse command line arguments
argv <- commandArgs(trailingOnly = TRUE)

SAMPLE_NAME <- argv[1]

cat("Processing cellsnp results for:", SAMPLE_NAME, "\n")

# Define input files
tsv_file <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a/cellSNP.samples.tsv")
AD_mtx <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a/cellSNP.tag.AD.mtx")
DP_mtx <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a/cellSNP.tag.DP.mtx")
OTH_mtx <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a/cellSNP.tag.OTH.mtx")
VCF_file <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a/cellSNP.base.vcf")
output_file <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a.summary.tsv.gz")

# Check if output already exists
if (file.exists(output_file)) {
    cat("Output file already exists:", output_file, "\n")
    quit(save = "no", status = 0)
}

# Check if all input files exist
input_files <- c(AD_mtx, DP_mtx, OTH_mtx, VCF_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
    cat("ERROR: Missing input files:\n")
    cat(paste(missing_files, collapse = "\n"), "\n")
    quit(save = "no", status = 1)
}

cat("Loading matrices...\n")

# Load matrices
AD <- Matrix::readMM(AD_mtx)
DP <- Matrix::readMM(DP_mtx)
OTH <- Matrix::readMM(OTH_mtx)
RD <- DP - AD  # Reference depth = Total depth - Alt depth

# Load VCF
cat("Loading VCF...\n")
vcf <- read_tsv(VCF_file, comment = "#", col_names = FALSE, show_col_types = FALSE)

# Load barcodes
barcode_file <- paste0("barcodes/", SAMPLE_NAME, "_barcodes.txt")
if (!file.exists(barcode_file)) {
    stop("Barcode file not found: ", barcode_file)
}
barcodes <- read.table(barcode_file, stringsAsFactors = FALSE)

cat("Processing", nrow(vcf), "SNPs and", nrow(barcodes), "barcodes...\n")

# Create SNP information
SNP_INFO <- vcf %>%
    mutate(SNP_ID = str_c(X1, ":", X2, ":", X4, ":", X5)) %>%
    dplyr::select(SNP_ID, X1, X2, X4, X5)
colnames(SNP_INFO) <- c("SNP_ID", "CHR", "POS", "REF", "ALT")

# Process matrices in chunks to handle memory efficiently
process_matrix_chunks <- function(matrix, matrix_name) {
    cat("Processing", matrix_name, "matrix...\n")
    
    # Convert to tibble and add column names
    result <- matrix %>% 
        as.matrix() %>% 
        as_tibble()
    colnames(result) <- barcodes$V1
    result$SNP <- SNP_INFO$SNP_ID
    
    # Process in chunks to reduce memory usage
    chunk_size <- 10000
    snp_chunks <- split(result, ceiling(seq_len(nrow(result)) / chunk_size))
    
    result_long <- lapply(snp_chunks, function(chunk) {
        # Select only columns with non-zero values
        chunk_filtered <- chunk %>%
            select(SNP, where(~ any(. > 0)))
        
        if (ncol(chunk_filtered) > 1) {
            chunk_filtered %>%
                pivot_longer(names_to = "cell_barcode", values_to = matrix_name, -SNP) %>%
                filter(.data[[matrix_name]] > 0)
        } else {
            NULL
        }
    }) %>%
        bind_rows()
    
    return(result_long)
}

# Process each matrix
ref_res_long <- process_matrix_chunks(RD, "REFcount")
alt_res_long <- process_matrix_chunks(AD, "ALTcount")
oth_res_long <- process_matrix_chunks(OTH, "OTHcount")

cat("Merging results...\n")

# Combine all results
result <- full_join(ref_res_long, alt_res_long, by = c("SNP", "cell_barcode")) %>%
    full_join(oth_res_long, by = c("SNP", "cell_barcode")) %>%
    left_join(SNP_INFO, by = c("SNP" = "SNP_ID"))

# Replace NA values with 0
result[is.na(result)] <- 0

# Remove duplicates
result <- distinct(result)

cat("Final result dimensions:", nrow(result), "x", ncol(result), "\n")

# Write output
cat("Writing output to:", output_file, "\n")
write_tsv(result, output_file)

cat("Processing completed successfully!\n")
