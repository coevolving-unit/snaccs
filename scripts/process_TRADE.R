library(TRADE)
library(dplyr)

# Read samples list and extract regions
samples_list <- read.csv('/data/region_samples_list.csv')
regions <- levels(as.factor(samples_list$region))

# Define cell type categories
subclass_other <- c('Oligo', 'OPC', 'Astro', 'Fibro', 'Micro', 'BEC', 'SM', 'Peri', 'Tcell')
subclass_inn <- c('ADARB2', 'Chandelier', 'VIP', 'SST', 'LAMP5', 'PVALB', 'PAX6')
subclass_exn <- c('L23IT', 'L4IT', 'L5IT', 'L6IT', 'L5ET', 'L56NP', 'L6b', 'L6CT')

process_de_results <- function(regions, cell_types, de_type, file_prefix, suffix) {
  res <- data.frame()
  for (region in regions) {
    print(region)
    for (cell_type in cell_types) {
      print(cell_type)
      file_path <- paste(region, '_', cell_type, '_', de_type, '_outliers.rds', sep = "")
      if (file.exists(file_path)) {
        loaded <- readRDS(file_path)
      		sig=data.frame(twi = loaded$distribution_summary$transcriptome_wide_impact,
      					region <- region,
      					celltype <- cell_type)
          res <- rbind(res, sig)
      } else {
        print(paste('Not analyzed:', file_path))
      }
    }
  }
  print(paste(file_prefix, '_', suffix,'_outliers.rds', sep = ""))
  saveRDS(res, file = paste(file_prefix, '_', suffix,'_outliers.rds', sep = ""))
}

# Process data for subclass
process_de_results(regions, subclass_other, 'TRADE_sex', 'other_SEX_TRADE','class')
process_de_results(regions, subclass_inn, 'TRADE_sex', 'inn_SEX_TRADE','class')
process_de_results(regions, subclass_exn, 'TRADE_sex', 'exn_SEX_TRADE','class')

process_de_results(regions, subclass_other, 'TRADE_age', 'other_age_TRADE','class')
process_de_results(regions, subclass_inn, 'TRADE_age', 'inn_age_TRADE','class')
process_de_results(regions, subclass_exn, 'TRADE_age', 'exn_age_TRADE','class')

## load and combine results

# sex + class
other_res = readRDS('other_SEX_TRADE_class_outliers.rds')
inn_res = readRDS('inn_SEX_TRADE_class_outliers.rds')
exn_res = readRDS('exn_SEX_TRADE_class_outliers.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_SEX_TRADE_class_outliers.rds') # done

# age + class
other_res = readRDS('other_age_TRADE_class_outliers.rds')
inn_res = readRDS('inn_age_TRADE_class_outliers.rds')
exn_res = readRDS('exn_age_TRADE_class_outliers.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_AGE_TRADE_class_outliers.rds') # done

# General function to process data
process_trade_data <- function(regions, cell_types, file_prefix, suffix) {
  res <- data.frame()
  for (region in regions) {
    print(region)
    for (cell_type in cell_types) {
      print(cell_type)
      file_path <- paste(region, cell_type, "TRADE-jk-pval.rds", sep = "_")
      if (file.exists(file_path)) {
        loaded <- readRDS(file_path)
        loaded = data.frame(p = loaded)
        loaded$region <- region
        loaded$celltype <- cell_type
        res <- rbind(res, loaded)
      } else {
        print(paste('Not analyzed:', file_path))
      }
    }
  }
  print(paste(file_prefix, '_','TRADE_jk_pval','_', suffix, '.rds', sep = ""))
  saveRDS(res, file = paste(file_prefix, '_','TRADE_jk_pval','_', suffix, '.rds', sep = ""))
}

# Process data for subclass
process_trade_data(regions, subclass_other, 'other', 'class') 
process_trade_data(regions, subclass_inn, 'inn', 'class') 
process_trade_data(regions, subclass_exn, 'exn', 'class') 

# General function to process data
process_trade_data <- function(regions, cell_types, file_prefix, suffix) {
  res <- data.frame()
  for (region in regions) {
    print(region)
    for (cell_type in cell_types) {
      print(cell_type)
      file_path <- paste(region, cell_type, "DE-loo-TRADE-jk.rds", sep = "_")
      if (file.exists(file_path)) {
        loaded <- readRDS(file_path)
        loaded = data.frame(loaded)
        loaded$region <- region
        loaded$celltype <- cell_type
        res <- rbind(res, loaded)
      } else {
        print(paste('Not analyzed:', file_path))
      }
    }
  }
  print(paste(file_prefix, '_','TRADE_loo','_', suffix, '.rds', sep = ""))
  saveRDS(res, file = paste(file_prefix, '_','TRADE_loo','_', suffix, '.rds', sep = ""))
}

# Process data for subclass
process_trade_data(regions, subclass_other, 'other', 'class') 
process_trade_data(regions, subclass_inn, 'inn', 'class') 
process_trade_data(regions, subclass_exn, 'exn', 'class') 

# sex + class
other_res = readRDS('other_TRADE_jk_pval_class.rds')
inn_res = readRDS('inn_TRADE_jk_pval_class.rds')
exn_res = readRDS('exn_TRADE_jk_pval_class.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_TRADE_jk_pval_class.rds') 

other_res = readRDS('other_TRADE_loo_class.rds')
inn_res = readRDS('inn_TRADE_loo_class.rds')
exn_res = readRDS('exn_TRADE_loo_class.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_TRADE_loo_class.rds')



