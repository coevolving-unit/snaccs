library(dplyr)

# Read samples list and extract regions
samples_list <- read.csv('/data/region_samples_list.csv')
regions <- levels(as.factor(samples_list$region))

# Define cell type categories
subclass_other <- c('Oligo', 'OPC', 'Astro', 'Fibro', 'Micro', 'BEC', 'SM', 'Peri', 'Tcell')
subclass_inn <- c('ADARB2', 'Chandelier', 'VIP', 'SST', 'LAMP5', 'PVALB', 'PAX6')
subclass_exn <- c('L23IT', 'L4IT', 'L5IT', 'L6IT', 'L5ET', 'L56NP', 'L6b', 'L6CT')

# General function to process data
process_data <- function(regions, cell_types, file_prefix, suffix) {
  res <- data.frame()
  for (region in regions) {
    print(region)
    for (cell_type in cell_types) {
      print(cell_type)
      file_path <- paste(region, cell_type, 'varPart_continuous_corrected.rds', sep = "_")
      if (file.exists(file_path)) {
        loaded <- readRDS(file_path)
        loaded$region <- region
        loaded$celltype <- cell_type
        loaded$gene <- rownames(loaded)
        rownames(loaded) <- NULL
        res <- rbind(res, loaded)
      } else {
        print(paste('Not analyzed:', file_path))
      }
    }
  }
  print(paste(file_prefix, '_','varpart_continuous_corrected', '_',suffix, '.rds', sep = ""))
  saveRDS(res, file = paste(file_prefix, '_','varpart_continuous_corrected', '_',suffix, '.rds', sep = ""))
}

# Process data for subclass
process_data(regions, subclass_other, 'other', 'class') 
process_data(regions, subclass_inn, 'inn', 'class') 
process_data(regions, subclass_exn, 'exn', 'class') 

other_res = readRDS('other_varpart_continuous_corrected_class.rds')
inn_res = readRDS('inn_varpart_continuous_corrected_class.rds')
exn_res = readRDS('exn_varpart_continuous_corrected_class.rds')
all_res = rbind(other_res, inn_res, exn_res)
saveRDS(all_res, file = 'all_varpart_continuous_corrected_class.rds') 

