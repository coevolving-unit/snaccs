print("loading data")

library(variancePartition)

arg = commandArgs(trailingOnly=TRUE)
region_now = arg[1]
celltype_now = arg[2]

ynow = readRDS(paste(region_now,celltype_now,"normalized_expression_corrected_outliers.rds",sep="_"))

meta = readRDS('meta_merged.rds')

c = meta[,c('new.id',celltype_now)]
colnames(c)[2] = 'count'
keep = subset(c, count >= 30)$new.id
ynow = ynow[,which(colnames(ynow) %in% keep)]

meta_data_now = subset(meta, new.id %in% colnames(ynow))
meta_data_now = unique(meta_data_now[,c('new.id','region','RIN','PMI','individual',
										'batch','sex','age','nreads','maprate','total_cells')])
rownames(meta_data_now) = meta_data_now$new.id
idx = match(colnames(ynow), meta_data_now$new.id)
meta_data_now = meta_data_now[idx,]
print(table(meta_data_now$sex))

print('scale all continuous variables')

meta_data_now$RIN = scale(meta_data_now$RIN)
meta_data_now$PMI = scale(meta_data_now$PMI)
meta_data_now$nreads = scale(meta_data_now$nreads)
meta_data_now$total_cells = scale(meta_data_now$total_cells)
meta_data_now$age = scale(meta_data_now$age)

print('variance partitioning')

design = as.formula(paste('~',paste(c('age','RIN','PMI','nreads','total_cells','(1|sex)'),collapse=' + ')))
varPart = fitExtractVarPartModel(ynow, design, meta_data_now)

saveRDS(varPart, file = paste(region_now,celltype_now,"varPart_corrected.rds",sep="_"))

meta_data_now$sex = ifelse(meta_data_now$sex == 'Male', 1, -1)

design = as.formula(paste('~',paste(c('age','RIN','PMI','nreads','total_cells','sex'),collapse=' + ')))
varPart = fitExtractVarPartModel(ynow, design, meta_data_now)

saveRDS(varPart, file = paste(region_now,celltype_now,"varPart_continuous_corrected.rds",sep="_"))

print('done')
