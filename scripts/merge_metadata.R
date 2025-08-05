#!/usr/bin/env Rscript

# scripts/merge_metadata.R
# Merge sample metadata from multiple sources

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
})

cat("Merging sample metadata from multiple sources...\n")

# Set up paths - customize these for your data
METADATA_DIR <- "metadata"
OUTPUTS_DIR <- "outputs"

# Ensure output directory exists
dir.create(OUTPUTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Load brain-level and sample-level metadata
# Contains PMI, pH, RIN, sample information
metadata_file <- file.path(METADATA_DIR, "region_samples_list.txt")
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file, 
       "\nPlease ensure region_samples_list.txt is in the metadata/ directory")
}

cat("Loading sample metadata...\n")
meta_data <- read.delim(metadata_file)

# Clean and process metadata
meta_data$RIN <- as.numeric(meta_data$RIN)
meta_data$orig.ident <- str_sub(meta_data$sample_id, -5, -1)
meta_data$new.id <- str_sub(meta_data$new_id, -5, -1)

# Load Cell Ranger sequencing metrics
seq_metrics_file <- file.path(METADATA_DIR, "sequencing_metrics.csv")
if (file.exists(seq_metrics_file)) {
  cat("Loading Cell Ranger sequencing metrics...\n")
  seq_info <- read.csv(seq_metrics_file)
  seq_info <- subset(seq_info, type == 'all')
  
  # Merge with main metadata
  meta_data <- merge(meta_data, seq_info, by = 'new.id')
  colnames(meta_data)[2] <- 'orig.ident'
  meta_data$orig.ident <- str_sub(meta_data$orig.ident, -5, -1)
  
  # Clean numeric columns
  meta_data$Number.of.Reads <- gsub(",", "", meta_data$Number.of.Reads)
  meta_data$Number.of.Reads <- as.numeric(meta_data$Number.of.Reads)
  
  meta_data$Reads.Mapped.to.Genome <- gsub("%", "", meta_data$Reads.Mapped.to.Genome)
  meta_data$Reads.Mapped.to.Genome <- as.numeric(meta_data$Reads.Mapped.to.Genome)
  
  meta_data$Estimated.Number.of.Cells <- gsub(",", "", meta_data$Estimated.Number.of.Cells)
  meta_data$Estimated.Number.of.Cells <- as.numeric(meta_data$Estimated.Number.of.Cells)
  
  cat("Cell Ranger metrics merged successfully\n")
} else {
  cat("Warning: Cell Ranger sequencing metrics not found at:", seq_metrics_file, "\n")
}

# Load Alevin-fry QC metrics
alevin_qc_file <- file.path(METADATA_DIR, "alevin_qc_tables.csv")
if (file.exists(alevin_qc_file)) {
  cat("Loading Alevin-fry QC metrics...\n")
  seq_info_alevin <- read.csv(alevin_qc_file)
  
  # Create summary dataframe
  alevin_summary <- data.frame(
    orig.ident = seq_info_alevin$sample,
    nreads = seq_info_alevin$Total.number.of.processed.reads,
    maprate = seq_info_alevin$Number.of.mapped.reads / seq_info_alevin$Total.number.of.processed.reads
  )
  alevin_summary$orig.ident <- str_sub(alevin_summary$orig.ident, -5, -1)
  
  cat("Alevin-fry QC metrics processed successfully\n")
} else {
  cat("Warning: Alevin-fry QC metrics not found at:", alevin_qc_file, "\n")
  # Create empty dataframe if file doesn't exist
  alevin_summary <- data.frame(orig.ident = character(0), nreads = numeric(0), maprate = numeric(0))
}

# Load cell count data from annotation results
cat("Loading cell count data...\n")
cell_count_files <- c(
  "other_annotate_class_counts_corrected.csv",
  "inn_annotate_class_counts_corrected.csv", 
  "exn_upper_annotate_class_counts_corrected.csv",
  "exn_lower_annotate_class_counts_corrected.csv"
)

# Check if count files exist
count_data_list <- list()
missing_files <- character(0)

for (file in cell_count_files) {
  filepath <- file.path(OUTPUTS_DIR, file)
  if (file.exists(filepath)) {
    count_data_list[[file]] <- read.csv(filepath)
    cat("Loaded:", file, "\n")
  } else {
    missing_files <- c(missing_files, file)
  }
}

# Combine cell count data if files exist
if (length(count_data_list) > 0) {
  # Combine all count data
  all_counts <- count_data_list[[1]]
  
  if (length(count_data_list) > 1) {
    for (i in 2:length(count_data_list)) {
      all_counts <- cbind(all_counts, count_data_list[[i]][, -1])
    }
  }
  
  # Calculate total cells
  count_cols <- 2:ncol(all_counts)
  all_counts$total_cells <- rowSums(all_counts[, count_cols])
  colnames(all_counts)[1] <- 'new.id'
  
  cat("Cell count data combined successfully\n")
} else {
  cat("Warning: No cell count files found. Creating empty dataframe.\n")
  all_counts <- data.frame(new.id = character(0), total_cells = numeric(0))
}

if (length(missing_files) > 0) {
  cat("Missing cell count files:\n")
  for (file in missing_files) {
    cat("  -", file, "\n")
  }
}

# Merge all metadata sources
cat("Merging all metadata sources...\n")

# Remove duplicate column if it exists
if ("sample_id" %in% colnames(meta_data)) {
  meta_data <- meta_data[, !colnames(meta_data) %in% "sample_id"]
}

# Merge datasets
meta_all <- meta_data
if (nrow(alevin_summary) > 0) {
  meta_all <- merge(meta_all, alevin_summary, by = 'orig.ident', all = TRUE)
}
if (nrow(all_counts) > 0) {
  meta_all <- merge(meta_all, all_counts, by = 'new.id', all = TRUE)
}

cat("Current metadata dimensions:", nrow(meta_all), "x", ncol(meta_all), "\n")

# Handle duplicate samples (replicates)
# For samples with multiple technical replicates, average mapping rates and sum read counts
if ("nreads" %in% colnames(meta_all) && "maprate" %in% colnames(meta_all)) {  
  cat("Processing technical replicates...\n")
  
  # Identify samples with replicates
  replicate_summary <- meta_all %>% 
    group_by(new.id) %>% 
    summarise(n = n(), .groups = 'drop')
  
  duplicated_samples <- subset(replicate_summary, n > 1)
  
  if (nrow(duplicated_samples) > 0) {
    cat("Found", nrow(duplicated_samples), "samples with replicates\n")
    
    # Calculate aggregated metrics for duplicated samples
    replicate_stats <- meta_all %>%
      filter(new.id %in% duplicated_samples$new.id) %>%
      group_by(new.id) %>%
      summarise(
        nreads = sum(nreads, na.rm = TRUE),
        maprate = mean(maprate, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Update metadata with aggregated stats
    for (i in 1:nrow(replicate_stats)) {
      sample_id <- replicate_stats$new.id[i]
      meta_all[meta_all$new.id == sample_id, c('nreads', 'maprate')] <- 
        replicate_stats[i, c('nreads', 'maprate')]
    }
    
    cat("Technical replicates processed successfully\n")
  }
}

# Handle specific sample corrections (customize as needed)
# Example: Fix batch assignment for specific samples
if ("batch" %in% colnames(meta_all)) {
  # Add specific sample corrections here
  meta_all[which(meta_all$orig.ident == 79557),'batch'] = meta_all[which(meta_all$orig.ident == 81708),'batch']
  cat("Sample-specific corrections applied\n")
}

# Save merged metadata
output_file <- file.path(OUTPUTS_DIR, "meta_merged.rds")
saveRDS(meta_all, file = output_file)

cat("Merged metadata saved to:", output_file, "\n")
cat("Final metadata dimensions:", nrow(meta_all), "x", ncol(meta_all), "\n")

# Print summary statistics
cat("\nMetadata Summary:\n")
cat("================\n")
if ("sex" %in% colnames(meta_all)) {
  cat("Sex distribution:\n")
  print(table(meta_all$sex, useNA = "ifany"))
}
if ("batch" %in% colnames(meta_all)) {
  cat("Batch distribution:\n")
  print(table(meta_all$batch, useNA = "ifany"))
}
if ("total_cells" %in% colnames(meta_all)) {
  cat("Total cells - Mean:", round(mean(meta_all$total_cells, na.rm = TRUE)), 
      "Median:", round(median(meta_all$total_cells, na.rm = TRUE)), "\n")
}
if ("RIN" %in% colnames(meta_all)) {
  cat("RIN - Mean:", round(mean(meta_all$RIN, na.rm = TRUE), 2), 
      "Median:", round(median(meta_all$RIN, na.rm = TRUE), 2), "\n")
}

cat("\nMetadata merging completed successfully!\n")

