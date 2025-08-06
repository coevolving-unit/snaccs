#!/usr/bin/env Rscript

# QC analysis of X chromosome genotyping results

library(tidyverse)
library(viridis)
library(patchwork)

`%!in%` <- Negate(`%in%`)

# Parse command line arguments
argv <- commandArgs(trailingOnly = TRUE)

SAMPLE_NAME <- argv[1]
cat("Running QC analysis for:", SAMPLE_NAME, "\n")

# Define input files
RNA_SNPCELL <- paste0("Xcalls/", SAMPLE_NAME, ".RNA.chrX.1a.summary.tsv.gz")
RNA_ANNOVAR <- paste0("annovar/", SAMPLE_NAME, ".RNA.chrX.snp.call.hg38_multianno.txt")

# Check if input files exist
if (!file.exists(RNA_SNPCELL)) {
    stop("SNP-cell data not found: ", RNA_SNPCELL)
}
if (!file.exists(RNA_ANNOVAR)) {
    stop("ANNOVAR annotation not found: ", RNA_ANNOVAR)
}

# Create output directories
dir.create("QC", showWarnings = FALSE)
dir.create("snps", showWarnings = FALSE)

cat("Loading data...\n")

# Load RNA SNP-cell data
rna <- read_tsv(RNA_SNPCELL, show_col_types = FALSE)

# Function to load ANNOVAR annotation
load_anno <- function(ANNOVAR) {
    anno <- read_tsv(ANNOVAR, show_col_types = FALSE) %>%
        mutate(SNP_ID = str_c(Chr, Start, Ref, Alt, sep = ":")) %>%
        dplyr::select(SNP_ID, Func.refGene, Gene.refGene, ALL.sites.2015_08)
    colnames(anno) <- c("SNP_ID", "Region", "Gene", "ALL_Freq")
    anno$ALL_Freq[is.na(anno$ALL_Freq)] <- 0
    return(anno)
}

# Load annotation
rna_anno <- load_anno(RNA_ANNOVAR)

# Function to merge annotation with SNP data
merge_anno <- function(data, anno) {
    df <- left_join(data, anno, by = "SNP_ID") %>%
        filter(!str_detect(Gene, pattern = ";")) %>%
        filter(!is.na(Gene)) %>%
        filter(Region %in% c("intronic", "UTR5", "UTR3", "exonic", "ncRNA_exonic", "ncRNA_intronic", "splicing"))
    return(df)
}

# Merge RNA data with annotation
rna_df <- merge_anno(rna, rna_anno)
cat("Genes found:", length(unique(rna_df$Gene)), "\n")

# Save all SNP data
write_tsv(rna_df, paste0("QC/", SAMPLE_NAME, "_RNA_QC_ALL_SNP_df.tsv.gz"))

# Filter out rare variants (frequency > 0 from 1000 Genomes Project)
rna_df <- filter(rna_df, ALL_Freq > 0.0)
cat("Genes after frequency filtering:", length(unique(rna_df$Gene)), "\n")

# Save filtered SNP data
write_tsv(rna_df, paste0("QC/", SAMPLE_NAME, "_RNA_QC_passed_SNP_df.tsv.gz"))

# Load XCI annotation
xci_file <- '/data/genotype/xci-balaton.csv'
if (file.exists(xci_file)) {
    xci <- read.csv(xci_file)
    colnames(xci)[1] <- 'Gene'
    rna_df <- merge(rna_df, xci, by = 'Gene', all.x = TRUE)
    cat("XCI annotation loaded\n")
} else {
    cat("WARNING: XCI annotation file not found:", xci_file, "\n")
    # Create dummy columns
    rna_df$Transcript.type <- NA
    rna_df$Domain.category <- NA
}

cat("Processing SNP-level data...\n")

# Process SNPs: group by cell, SNP, and gene
snps <- rna_df %>% 
    group_by(cell_barcode, SNP_ID, Gene, Transcript.type, Domain.category) %>% 
    summarise(
        ref = sum(REFcount), 
        alt = sum(ALTcount), 
        other = sum(OTHcount),
        .groups = 'drop'
    )

snps$tot <- snps$ref + snps$alt + snps$other
snps$br <- snps$ref / snps$tot

# Define biallelic expression: SNP is biallelic in cell if 0.1 < br < 0.9
snps$bi <- ifelse(snps$br > 0.9, 'mono', 'bi')
snps$bi <- ifelse(snps$br < 0.1, 'mono', snps$bi)

cat("SNP biallelic status:\n")
print(table(snps$bi))

# Remove SNPs with insufficient reads (< 2 total reads)
snps_sub <- subset(snps, tot >= 2)
cat("SNPs after read depth filtering:", nrow(snps_sub), "\n")

cat("Processing gene-level data...\n")

# Gene-level summary
genes <- snps_sub %>% 
    group_by(Gene, Transcript.type, Domain.category) %>% 
    summarise(
        tot_n_snps = n(), 
        bi_n_snps = sum(bi == "bi"), 
        tot_bi_reads = sum(tot[bi == "bi"]),
        .groups = 'drop'
    )

genes$br <- genes$bi_n_snps / genes$tot_n_snps
genes$any.evidence <- ifelse(genes$bi_n_snps == 0, 'no', 'yes')
genes <- genes[order(genes$br), ]

cat("Genes with biallelic evidence:\n")
print(subset(genes, any.evidence == 'yes'))

# Load cell type annotations
cat("Loading cell type annotations...\n")
load_cell_annotations <- function(sample_name) {
    # Load different cell type annotation files
    annotation_files <- list(
        other = '/data/other_annotate_meta_corrected.csv',
        inn = '/data/inn_annotate_meta_corrected.csv',
        exnu = '/data/exn_upper_annotate_meta_corrected.csv',
        exnl = '/data/exn_lower_annotate_meta_corrected.csv'
    )
    
    anno_list <- list()
    sample_short <- str_sub(sample_name, -5, -1)
    
    for (name in names(annotation_files)) {
        file_path <- annotation_files[[name]]
        if (file.exists(file_path)) {
            df <- read.csv(file_path)
            df_subset <- subset(df, str_sub(orig.ident, 1, 5) == sample_short)
            if (nrow(df_subset) > 0) {
                anno_list[[name]] <- data.frame(
                    cell_barcode = df_subset$cell_id,
                    Annotation = df_subset$cell_class
                )
            }
        }
    }
    
    if (length(anno_list) > 0) {
        anno <- do.call(rbind, anno_list)
        anno$cell_barcode <- str_sub(anno$cell_barcode, 1, 16)
        return(anno)
    } else {
        return(NULL)
    }
}

# Load cell annotations
anno <- load_cell_annotations(SAMPLE_NAME)

if (!is.null(anno)) {
    cat("Cell type annotations loaded:", nrow(anno), "cells\n")
    
    # Merge with SNP data
    snps_sub$cell_barcode <- str_sub(snps_sub$cell_barcode, 1, 16)
    snps_sub_anno <- merge(snps_sub, anno, by = 'cell_barcode')
    
    # Gene-level analysis by cell type
    genes_by_celltype <- snps_sub_anno %>% 
        group_by(Gene, Annotation, Transcript.type, Domain.category) %>% 
        summarise(
            tot_n_snps = n(), 
            bi_n_snps = sum(bi == "bi"), 
            tot_bi_reads = sum(tot[bi == "bi"]),
            .groups = 'drop'
        )
    
    genes_by_celltype$br <- genes_by_celltype$bi_n_snps / genes_by_celltype$tot_n_snps
    genes_by_celltype$any.evidence <- ifelse(genes_by_celltype$bi_n_snps == 0, 'no', 'yes')
    genes_by_celltype$new_id <- SAMPLE_NAME
    
    # Cell-level analysis by cell type
    cells_by_celltype <- snps_sub_anno %>%
        group_by(Gene, Annotation, Domain.category) %>%
        summarise(
            total_cells = n_distinct(cell_barcode),
            cells_with_bi = n_distinct(cell_barcode[bi == "bi"]),
            cells_only_mono = n_distinct(cell_barcode[bi == "mono" & 
                !(cell_barcode %in% cell_barcode[bi == "bi"])]),
            .groups = 'drop'
        )
    cells_by_celltype$new_id <- SAMPLE_NAME
    
    # Save results with cell type information
    saveRDS(genes_by_celltype, paste0("snps/", SAMPLE_NAME, "_snp_counts_per_gene_celltype.rds"))
    saveRDS(cells_by_celltype, paste0("snps/", SAMPLE_NAME, "_cell_counts_per_gene_celltype.rds"))
    
    # Filter out cells with high XIST biallelic expression (potential doublets)
    cat("Filtering potential doublets based on XIST expression...\n")
    xist_bi <- subset(snps_sub_anno, Gene == 'XIST' & bi == 'bi')
    
    if (nrow(xist_bi) > 0) {
        cat("Found", length(unique(xist_bi$cell_barcode)), "cells with biallelic XIST expression\n")
        
        # Filter out these cells and reanalyze
        snps_sub_anno_filtered <- subset(snps_sub_anno, !(cell_barcode %in% unique(xist_bi$cell_barcode)))
        
        # Re-run gene-level analysis
        genes_filtered <- snps_sub_anno_filtered %>% 
            group_by(Gene, Annotation, Transcript.type, Domain.category) %>% 
            summarise(
                tot_n_snps = n(), 
                bi_n_snps = sum(bi == "bi"), 
                tot_bi_reads = sum(tot[bi == "bi"]),
                .groups = 'drop'
            )
        genes_filtered$br <- genes_filtered$bi_n_snps / genes_filtered$tot_n_snps
        genes_filtered$any.evidence <- ifelse(genes_filtered$bi_n_snps == 0, 'no', 'yes')
        genes_filtered$new_id <- SAMPLE_NAME
        
        # Re-run cell-level analysis
        cells_filtered <- snps_sub_anno_filtered %>%
            group_by(Gene, Annotation, Domain.category) %>%
            summarise(
                total_cells = n_distinct(cell_barcode),
                cells_with_bi = n_distinct(cell_barcode[bi == "bi"]),
                cells_only_mono = n_distinct(cell_barcode[bi == "mono" & 
                    !(cell_barcode %in% cell_barcode[bi == "bi"])]),
                .groups = 'drop'
            )
        cells_filtered$new_id <- SAMPLE_NAME
        
        # Save filtered results
        saveRDS(genes_filtered, paste0("snps/", SAMPLE_NAME, "_snp_counts_per_gene_celltype_filtered.rds"))
        saveRDS(cells_filtered, paste0("snps/", SAMPLE_NAME, "_cell_counts_per_gene_celltype_filtered.rds"))
        
        cat("Filtered results saved\n")
    }
    
} else {
    cat("WARNING: No cell type annotations found for", SAMPLE_NAME, "\n")
    
    # Save basic results without cell type information
    genes$new_id <- SAMPLE_NAME
    saveRDS(genes, paste0("snps/", SAMPLE_NAME, "_snp_counts_per_gene.rds"))
}

cat("QC analysis completed for", SAMPLE_NAME, "\n")
cat("Results saved to:\n")
cat("  - QC/", SAMPLE_NAME, "_RNA_QC_ALL_SNP_df.tsv.gz\n")
cat("  - QC/", SAMPLE_NAME, "_RNA_QC_passed_SNP_df.tsv.gz\n")
cat("  - snps/", SAMPLE_NAME, "_*.rds\n")
