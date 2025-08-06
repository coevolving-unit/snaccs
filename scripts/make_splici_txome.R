#!/usr/bin/env Rscript

# scripts/make_splici_txome.R
# Creates splici transcriptomes for both male and female references

suppressPackageStartupMessages({
  library(eisaR)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
  library(GenomicFeatures)
})

cat("Creating splici transcriptomes...\n")

# Parameters
read_length <- 91
flank_trim_length <- 5

# Create output directories
dir.create("data/references/transcriptomes", recursive = TRUE, showWarnings = FALSE)

# Male genome (with Y chromosome)
cat("Creating male splici transcriptome...\n")
gtf_path_male <- "data/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
genome_path_male <- "data/references/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
output_dir_male <- paste0("data/references/transcriptomes/transcriptome_splici_fl", read_length - flank_trim_length)

if (file.exists(gtf_path_male) && file.exists(genome_path_male)) {
  make_splici_txome(gtf_path = gtf_path_male, 
                    genome_path = genome_path_male, 
                    read_length = read_length, 
                    flank_trim_length = flank_trim_length, 
                    output_dir = output_dir_male)
  cat("Male splici transcriptome created successfully at:", output_dir_male, "\n")
} else {
  cat("Male reference files not found:\n")
  cat("GTF:", gtf_path_male, "exists:", file.exists(gtf_path_male), "\n")
  cat("Genome:", genome_path_male, "exists:", file.exists(genome_path_male), "\n")
}

# Female genome (Y chromosome masked)
cat("Creating female splici transcriptome...\n")
gtf_path_female <- "data/references/GRCh38_noY/genes/genes.gtf"
genome_path_female <- "data/references/GRCh38_noY/fasta/genome.fa"
output_dir_female <- paste0("data/references/transcriptomes/transcriptome_splici_fl_noY", read_length - flank_trim_length)

if (file.exists(gtf_path_female) && file.exists(genome_path_female)) {
  make_splici_txome(gtf_path = gtf_path_female, 
                    genome_path = genome_path_female, 
                    read_length = read_length, 
                    flank_trim_length = flank_trim_length, 
                    output_dir = output_dir_female)
  cat("Female splici transcriptome created successfully at:", output_dir_female, "\n")
} else {
  cat("Female reference files not found:\n")
  cat("GTF:", gtf_path_female, "exists:", file.exists(gtf_path_female), "\n")
  cat("Genome:", genome_path_female, "exists:", file.exists(genome_path_female), "\n")
}

cat("Splici transcriptome creation completed\n")
cat("Output files saved to: data/references/transcriptomes/\n")
