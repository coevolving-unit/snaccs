# scripts/utils/load_fry.R
# Function to load Alevin-fry quantification results

load_fry <- function(frydir, which_counts = c('S', 'U', 'A'), verbose = FALSE) {
  suppressPackageStartupMessages({
    library(rjson)
    library(Matrix)
    library(SingleCellExperiment)
  })
  
  # Validate inputs
  if (!dir.exists(frydir)) {
    stop("Directory does not exist: ", frydir)
  }
  
  # Read metadata
  qfile <- file.path(frydir, "quant.json")
  if (!file.exists(qfile)) {
    qfile <- file.path(frydir, "meta_info.json")
  }
  
  if (!file.exists(qfile)) {
    stop("Neither quant.json nor meta_info.json found in ", frydir)
  }
  
  meta_info <- fromJSON(file = qfile)
  ng <- meta_info$num_genes
  usa_mode <- meta_info$usa_mode
  
  if (usa_mode) {
    if (length(which_counts) == 0) {
      stop("Please provide at least one status in 'U', 'S', 'A'")
    }
    valid_counts <- c('U', 'S', 'A')
    if (!all(which_counts %in% valid_counts)) {
      stop("Invalid count types. Must be one of: ", paste(valid_counts, collapse = ", "))
    }
    if (verbose) {
      message("Processing in USA mode, returning ", paste(which_counts, collapse = '+'))
    }
  } else if (verbose) {
    message("Processing in standard mode, returning spliced counts")
  }
  
  # Read count matrix
  mtx_file <- file.path(frydir, "alevin", "quants_mat.mtx")
  if (!file.exists(mtx_file)) {
    stop("Count matrix not found: ", mtx_file)
  }
  
  af_raw <- readMM(mtx_file)
  
  # Handle USA mode
  if (usa_mode) {
    if (ng %% 3 != 0) {
      stop("Number of quantified targets is not a multiple of 3")
    }
    ng <- as.integer(ng / 3)
  }
  
  # Read gene and cell files
  gene_file <- file.path(frydir, "alevin", "quants_mat_cols.txt")
  cell_file <- file.path(frydir, "alevin", "quants_mat_rows.txt")
  
  if (!file.exists(gene_file) || !file.exists(cell_file)) {
    stop("Gene or cell file missing")
  }
  
  afg <- read.csv(gene_file, strip.white = TRUE, header = FALSE, 
                  nrows = ng, col.names = "gene_ids", row.names = 1)
  afc <- read.csv(cell_file, strip.white = TRUE, header = FALSE,
                  col.names = "barcodes", row.names = 1)
  
  # Sum counts according to which_counts
  if (usa_mode) {
    rd <- list("S" = seq(1, ng), 
               "U" = seq(ng + 1, 2 * ng),
               "A" = seq(2 * ng + 1, 3 * ng))
    
    o <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
    if (length(which_counts) > 1) {
      for (wc in which_counts[-1]) {
        o <- o + af_raw[, rd[[wc]], drop = FALSE]
      }
    }
  } else {
    o <- af_raw
  }
  
  # Create SingleCellExperiment
  sce <- SingleCellExperiment(
    list(counts = t(o)),
    colData = afc,
    rowData = afg
  )
  
  return(sce)
}

---

#!/usr/bin/env Rscript

# scripts/alevin_matrix_M.R
# Create SingleCellExperiment objects for male samples

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
fry_path <- file.path("mapped_reads_M", quant_res)

# Load data
af <- load_fry(fry_path)

# Save SingleCellExperiment object
output_file <- file.path("mapped_reads_M", paste0(sample_id, "_sce.rds"))
saveRDS(af, file = output_file)

cat("Saved SCE object to:", output_file, "\n")
cat("Done processing sample:", sample_id, "\n")

---

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
