library(dplyr)

# Read samples list and extract regions
samples_list <- read.csv('/data/region_samples_list.csv')
regions <- levels(as.factor(samples_list$region))

# Define cell type categories
subclass_other <- c('Oligo', 'OPC', 'Astro', 'Micro', 'BEC')
subclass_inn <- c('ADARB2', 'Chandelier', 'VIP', 'SST', 'LAMP5', 'PVALB', 'PAX6')
subclass_exn <- c('L23IT', 'L4IT', 'L5IT', 'L6IT', 'L56NP', 'L6b', 'L6CT')

# Function to merge normalized expression data
merge_normalized_expression <- function(regions, cell_types, file_prefix, suffix) {
  res <- data.frame()
  for (region in regions) {
    for (cell_type in cell_types) {
      file_path <- paste(region, cell_type, 'normalized_expression_corrected_outliers.rds', sep = "_")
      if (file.exists(file_path)) {
        loaded <- readRDS(file_path)
        colnames(loaded) <- paste(colnames(loaded), region, cell_type, sep = ".")
        if (nrow(res) == 0) {
          res <- loaded
        } else {
          res <- merge(res, loaded, by = 'row.names', all = TRUE)
          rownames(res) <- res$Row.names
          res <- res[,-1]
        }
      } else {
        print(paste('Not analyzed:', file_path))
      }
    }
  }
  print(paste(file_prefix, 'normalized_expression',suffix,'corrected.rds', sep = "_"))
  saveRDS(res, file = paste(file_prefix, 'normalized_expression',suffix,'corrected_outliers.rds', sep = "_"))
}

# Function to process DE results
process_de_results <- function(regions, cell_types, de_type, file_prefix, suffix) {
  res <- data.frame()
  for (region in regions) {
    for (cell_type in cell_types) {
      file_path <- paste(region, cell_type, de_type, 'corrected.rds', sep = "_")
      if (file.exists(file_path)) {
        loaded <- readRDS(file_path)
        loaded$gene <- rownames(loaded)
        rownames(loaded) <- NULL
        loaded$region <- region
        loaded$celltype <- cell_type
        loaded$padj <- p.adjust(loaded$p, method = 'BH')
        res <- rbind(res, loaded)
      } else {
        print(paste('Not analyzed:', file_path))
      }
    }
  }
  print(paste(file_prefix, suffix, 'corrected.rds', sep = "_"))
  saveRDS(res, file = paste(file_prefix, suffix, 'corrected.rds', sep = "_"))
}

# Process data for subclass
merge_normalized_expression(regions, subclass_other, 'other', 'class') # done
merge_normalized_expression(regions, subclass_inn, 'inn', 'class') # done
merge_normalized_expression(regions, subclass_exn, 'exn', 'class') # done

process_de_results(regions, subclass_other, 'SEX_DE_limma_outliers', 'other_SEX_res_outliers', 'class') # done
process_de_results(regions, subclass_inn, 'SEX_DE_limma_outliers', 'inn_SEX_res_outliers', 'class') # done
process_de_results(regions, subclass_exn, 'SEX_DE_limma_outliers', 'exn_SEX_res_outliers', 'class') # done

process_de_results(regions, subclass_other, 'AGE_DE_limma_outliers', 'other_AGE_res_outliers', 'class') # done
process_de_results(regions, subclass_inn, 'AGE_DE_limma_outliers', 'inn_AGE_res_outliers', 'class') # done
process_de_results(regions, subclass_exn, 'AGE_DE_limma_outliers', 'exn_AGE_res_outliers', 'class') # done

# normalized expression per class
other_res = readRDS('other_normalized_expression_class_corrected_outliers.rds')
inn_res = readRDS('inn_normalized_expression_class_corrected_outliers.rds')
exn_res = readRDS('exn_normalized_expression_class_corrected_outliers.rds')
all_res = merge(other_res, inn_res, by = 'row.names', all = T)
rownames(all_res) = all_res$Row.names
all_res = all_res[,-1]
all_res = merge(all_res, exn_res, by = 'row.names', all = T)
rownames(all_res) = all_res$Row.names
all_res = all_res[,-1]
saveRDS(all_res, file = 'all_res_normalized_expression_class_corrected_outliers.rds') # done

# sex + class
other_res = readRDS('other_SEX_res_outliers_class_corrected.rds')
inn_res = readRDS('inn_SEX_res_outliers_class_corrected.rds')
exn_res = readRDS('exn_SEX_res_outliers_class_corrected.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_res_outliers_SEX_class_corrected.rds') # done

# age + class
other_res = readRDS('other_AGE_res_outliers_class_corrected.rds')
inn_res = readRDS('inn_AGE_res_outliers_class_corrected.rds')
exn_res = readRDS('exn_AGE_res_outliers_class_corrected.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_res_outliers_AGE_class_corrected.rds') # done

