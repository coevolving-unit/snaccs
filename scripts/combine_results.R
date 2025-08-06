#!/usr/bin/env Rscript

# Combine QC results from all samples

library(stringr)
library(dplyr)
library(tidyr)

cat("Combining QC results from all samples...\n")

# Load sample metadata
meta_file <- '/data/meta_merged.rds'
if (file.exists(meta_file)) {
    meta <- readRDS(meta_file)
    cat("Loaded metadata for", nrow(meta), "samples\n")
} else {
    cat("WARNING: Metadata file not found:", meta_file, "\n")
    meta <- NULL
}

# Function to combine gene-level results
combine_gene_results <- function(pattern, output_file) {
    cat("Combining files with pattern:", pattern, "\n")
    
    # Find all matching files
    samples <- list.files("snps/", pattern = pattern, full.names = TRUE)
    cat("Found", length(samples), "files\n")
    
    if (length(samples) == 0) {
        cat("No files found matching pattern:", pattern, "\n")
        return(NULL)
    }
    
    # Combine all data
    combined_data <- data.frame()
    
    for (i in 1:length(samples)) {
        cat("Processing:", samples[i], "\n")
        
        tryCatch({
            sample_data <- readRDS(samples[i])
            
            # Extract sample ID and add metadata
            if ('new_id' %in% colnames(sample_data)) {
                sample_data$orig.ident <- str_sub(sample_data$new_id, -5, -1)
            } else {
                # Extract from filename
                sample_name <- basename(samples[i])
                sample_id <- str_extract(sample_name, "sample_[0-9]+")
                sample_data$orig.ident <- str_sub(sample_id, -5, -1)
            }
            
            # Add metadata if available
            if (!is.null(meta) && 'orig.ident' %in% colnames(meta)) {
                sample_data <- merge(sample_data, 
                                   meta[, c('orig.ident', 'individual', 'region', 'sex')], 
                                   by = 'orig.ident', all.x = TRUE)
            }
            
            combined_data <- rbind(combined_data, sample_data)
            
        }, error = function(e) {
            cat("ERROR processing", samples[i], ":", e$message, "\n")
        })
    }
    
    if (nrow(combined_data) > 0) {
        cat("Combined data dimensions:", nrow(combined_data), "x", ncol(combined_data), "\n")
        write.csv(combined_data, file = output_file, row.names = FALSE)
        cat("Results saved to:", output_file, "\n")
        return(combined_data)
    } else {
        cat("No data to combine\n")
        return(NULL)
    }
}

# Combine gene-level results
cat("\n=== COMBINING GENE-LEVEL RESULTS ===\n")
genedat <- combine_gene_results("snp_counts_per_gene_celltype_filtered\\.rds", 'genedat.csv')

# Combine cell-level results
cat("\n=== COMBINING CELL-LEVEL RESULTS ===\n")
celldat <- combine_gene_results("cell_counts_per_gene_celltype_filtered\\.rds", 'celldat.csv')

# Generate summary statistics
if (!is.null(genedat)) {
    cat("\n=== GENE-LEVEL ANALYSIS ===\n")
    
    # Filter genes with sufficient data
    genedat_filtered <- subset(genedat, tot_n_snps > 2)
    
    cat("Summary of biallelic evidence by sex:\n")
    if ('sex' %in% colnames(genedat_filtered)) {
        summary_table <- table(genedat_filtered$any.evidence, genedat_filtered$sex)
        print(summary_table)
        
        # Genes with biallelic evidence
        biallelic_genes <- subset(genedat_filtered, any.evidence == 'yes')
        cat("\nGenes with biallelic evidence:\n")
        print(biallelic_genes[, c('Gene', 'Annotation', 'Domain.category', 'bi_n_snps', 'sex')])
        
        # Summary by domain category
        if ('Domain.category' %in% colnames(genedat_filtered)) {
            cat("\nSummary by XCI domain category:\n")
            domain_summary <- genedat_filtered %>% 
                group_by(sex, Domain.category, any.evidence) %>% 
                summarise(n = n(), .groups = 'drop')
            print(domain_summary)
        }
        
        # Compare male vs female
        cat("\nComparing male vs female results...\n")
        
        if ('Annotation' %in% colnames(genedat_filtered)) {
            # Create comparison table
            male_data <- subset(genedat_filtered, sex == 'Male') %>%
                group_by(Gene, Annotation, Domain.category, any.evidence) %>%
                summarise(n_male = n(), .groups = 'drop')
            
            female_data <- subset(genedat_filtered, sex == 'Female') %>%
                group_by(Gene, Annotation, Domain.category, any.evidence) %>%
                summarise(n_female = n(), .groups = 'drop')
            
            comparison <- merge(male_data, female_data, 
                              by = c("Gene", "Annotation", "Domain.category", "any.evidence"), 
                              all = TRUE)
            comparison$n_male[is.na(comparison$n_male)] <- 0
            comparison$n_female[is.na(comparison$n_female)] <- 0
            comparison$diff_f_minus_m <- comparison$n_female - comparison$n_male
            
            write.csv(comparison, file = 'genedat_comparison.csv', row.names = FALSE)
            cat("Male vs female comparison saved to: genedat_comparison.csv\n")
            
            # Show genes with biallelic evidence
            biallelic_comparison <- subset(comparison, any.evidence == 'yes')
            if (nrow(biallelic_comparison) > 0) {
                cat("\nGenes with biallelic evidence - sample counts:\n")
                print(biallelic_comparison)
            }
        }
    }
}

# Analyze specific genes of interest
if (!is.null(genedat)) {
    cat("\n=== SPECIFIC GENE ANALYSIS ===\n")
    
    # Genes of interest 
    genes_of_interest <- c('XIST', 'TSIX', 'FRMPD4', 'TBL1X', 'DMD', 'PTCHD1-AS', 
                          'IDS', 'FTX', 'DANT2', 'IL1RAPL1')
    
    available_genes <- intersect(genes_of_interest, unique(genedat$Gene))
    cat("Analyzing", length(available_genes), "genes of interest:", 
        paste(available_genes, collapse = ", "), "\n")
    
    if (length(available_genes) > 0) {
        gene_analysis <- subset(genedat, Gene %in% available_genes)
        
        if ('sex' %in% colnames(gene_analysis)) {
            # Summary by gene and sex
            gene_summary <- gene_analysis %>%
                group_by(Gene, sex, any.evidence) %>%
                summarise(n_samples = n(), 
                         mean_bi_snps = mean(bi_n_snps, na.rm = TRUE),
                         .groups = 'drop')
            
            write.csv(gene_summary, file = 'genes_of_interest_summary.csv', row.names = FALSE)
            cat("Gene-specific analysis saved to: genes_of_interest_summary.csv\n")
            
            # Show genes with sex differences
            cat("\nGenes with biallelic evidence by sex:\n")
            biallelic_by_sex <- subset(gene_summary, any.evidence == 'yes')
            print(biallelic_by_sex)
        }
    }
}

cat("\n=== ANALYSIS COMPLETED ===\n")
cat("Output files generated:\n")
cat("  - genedat.csv: Combined gene-level results\n")
cat("  - celldat.csv: Combined cell-level results\n")
if (file.exists('genedat_comparison.csv')) {
    cat("  - genedat_comparison.csv: Male vs female comparison\n")
}
if (file.exists('genes_of_interest_summary.csv')) {
    cat("  - genes_of_interest_summary.csv: Analysis of key XCI genes\n")
}
