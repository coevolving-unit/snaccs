library(variancePartition)
`%!in%` = Negate(`%in%`)

# load data

meta = readRDS('meta_merged.rds')

pseudo = readRDS('all_res_normalized_expression_class_corrected_outliers.rds') 

cells = sub('.*\\.', '', colnames(pseudo))
samples = sub("\\..*", "", colnames(pseudo))
samplecells = data.frame(new.id = samples, celltype = cells)
meta_data_now = subset(meta, new.id %in% samples)
meta_data_now = unique(meta_data_now[,c('new.id','region','RIN','PMI','individual',
										'batch','sex','age','nreads','maprate','total_cells')])
meta_data2 = merge(meta_data_now, samplecells, by = 'new.id', all.x = T)
meta_data2$key = paste(meta_data2$new.id, meta_data2$region,meta_data2$celltype, sep = ".")

check = reshape2::melt(meta[,c(3,12,41:151)])
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
meta_data2$key2 = paste(meta_data2$region,meta_data2$celltype)
meta_data2 = subset(meta_data2, key2 %!in% remove$key)

# match metadata

pseudo = pseudo[,which(colnames(pseudo) %in% meta_data2$key)]
idx = match(colnames(pseudo), meta_data2$key)
meta_data2 = meta_data2[idx,]

saveRDS(meta_data2, 'meta_all_res_normalized_expression_class_outliers.rds')

print('scale all continuous variables')

meta_data2$RIN = scale(meta_data2$RIN)
meta_data2$PMI = scale(meta_data2$PMI)
meta_data2$nreads = scale(meta_data2$nreads)
meta_data2$total_cells = scale(meta_data2$total_cells)
meta_data2$age = scale(meta_data2$age)
meta_data2$individual = as.character(meta_data2$individual)
rownames(meta_data2) = meta_data2$key

print('variance partitioning')

# limit to genes expressed in all regions x cell types
pseudo2 = pseudo[complete.cases(pseudo),]
dim(pseudo2) 

design = as.formula(paste('~',paste(c('age','RIN','PMI','nreads','total_cells',
									'(1|sex)','(1|region)','(1|celltype)','(1|individual)'),
									collapse=' + ')))
varPart = fitExtractVarPartModel(pseudo2, design, meta_data2)

saveRDS(varPart, file = "varpart_all_samples_class_corrected_outliers.rds") 
