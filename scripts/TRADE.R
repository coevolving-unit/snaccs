library(mashr)
library(TRADE)

print("loading data")

arg = commandArgs(trailingOnly=TRUE)
region_now = arg[1]
celltype_now = arg[2]

loaded = readRDS(paste(region_now,celltype_now,'SEX_DE_limma_outliers_corrected.rds',sep = "_"))
loaded$se = sqrt(loaded$s2) * loaded$stdev
colnames(loaded) = c('log2FoldChange','s2','stdev','pvalue','lfcSE')
TRADE_output = TRADE(mode = "univariate",
                     results1 = loaded,
                     annot_table = NULL,
                     genes_exclude = NULL,
                     n_sample = NULL)
                     
saveRDS(TRADE_output, paste(region_now,celltype_now,'TRADE_sex_outliers.rds',sep = "_"))
                     
loaded = readRDS(paste(region_now,celltype_now,'AGE_DE_limma_outliers_corrected.rds',sep = "_"))
loaded$se = sqrt(loaded$s2) * loaded$stdev
colnames(loaded) = c('log2FoldChange','s2','stdev','pvalue','lfcSE')
TRADE_output = TRADE(mode = "univariate",
                     results1 = loaded,
                     annot_table = NULL,
                     genes_exclude = NULL,
                     n_sample = NULL)
                     
saveRDS(TRADE_output, paste(region_now,celltype_now,'TRADE_age_outliers.rds',sep = "_"))

print('done')        
