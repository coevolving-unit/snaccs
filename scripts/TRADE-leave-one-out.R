library(Seurat)
library(SeuratObject)
library(WGCNA)
library(Matrix)
library(DescTools)
library(glmGamPoi)
library(slanter)
library(harmony)
library(scRNAseq)
library(ExperimentHub)
library(dplyr)
library(limma)
library(edgeR)
library(SingleCellExperiment)
library(Matrix.utils)
library(stringr)
library(mashr)
library(TRADE)
library(DESeq2)

`%!in%` = Negate(`%in%`)

print("loading data")

arg = commandArgs(trailingOnly=TRUE)
region_now = arg[1]
celltype_now = arg[2]

pseudo = readRDS(paste(region_now,celltype_now,'pseudo_corrected.rds',sep="-"))

meta = readRDS('meta_merged.rds')
colnames(meta)[which(colnames(meta) == 'BEC.arterial')] = 'arterial'
colnames(meta)[which(colnames(meta) == 'BEC.venous')] = 'venous'
colnames(meta)[which(colnames(meta) == 'BEC.cap')] = 'cap'

# remove samples with < 30 cells in cell type

c = meta[,c('new.id',celltype_now)]
colnames(c)[2] = 'count'
keep = subset(c, count >= 30)$new.id

pseudo = pseudo[,which(colnames(pseudo) %in% keep)]

meta_data_now = subset(meta, new.id %in% colnames(pseudo))
meta_data_now = unique(meta_data_now[,c('new.id','region','RIN','PMI','ph','individual',
										'batch','sex','age','nreads','maprate','total_cells')])
idx = match(colnames(pseudo), meta_data_now$new.id)
meta_data_now = meta_data_now[idx,]
print(table(meta_data_now$sex))
	
# remove mito and ribosomal genes

rpl = rownames(pseudo)[which(str_sub(rownames(pseudo), 1, 3) == 'RPL')]
rps = rownames(pseudo)[which(str_sub(rownames(pseudo), 1, 3) == 'RPS')]
mrpl = rownames(pseudo)[which(str_sub(rownames(pseudo), 1, 4) == 'MRPL')]
mrps = rownames(pseudo)[which(str_sub(rownames(pseudo), 1, 4) == 'MRPS')]
ribo = c(rpl, rps, mrpl, mrps)
ribo = ribo[which(ribo %!in% c('RPS4X', 'RPS4Y1', 'RPS4Y2'))]

mt = rownames(pseudo)[which(str_sub(rownames(pseudo), 1, 3) == 'MT-')]

pseudo = pseudo[which(rownames(pseudo) %!in% ribo),]
pseudo = pseudo[which(rownames(pseudo) %!in% mt),]

# remove lowly expressed genes

cpm_cutoff = 5
cpm_means = rowMeans(cpm(pseudo))
keep = names(cpm_means[which(cpm_means > cpm_cutoff)])
print(length(keep))

# subsample

out = data.frame()
for(i in 1:length(meta_data_now$new.id)){

	print(paste(i," out of ",length(meta_data_now$new.id),sep=""))
	snow = meta_data_now$new.id[i]
	meta_data_now_sub = meta_data_now[which(meta_data_now$new.id != snow),]
	pseudo_sub = pseudo[,which(colnames(pseudo) != snow)]

	f = length(subset(meta_data_now_sub, sex == 'Female')$sex)
	m = length(subset(meta_data_now_sub, sex == 'Male')$sex)
	
	print('normalization')
	
	d0 = DGEList(pseudo_sub)
	d0 = d0[keep, , keep.lib.sizes=FALSE]
	d0 = calcNormFactors(d0, method = 'TMM')
	dds = DESeqDataSetFromMatrix(countData = d0$counts,
                              colData = meta_data_now_sub,
                              design = ~ sex)
	dds = DESeq(dds)
	cooks = assays(dds)[["cooks"]]
	cooks_outlier = cooks > qf(0.99, df1 = 2, df2 = ncol(dds) - 2)
	outlier_genes = rownames(dds)[rowSums(cooks_outlier) > 0]
	dds_clean = replaceOutliers(dds)
	d0$counts = counts(dds_clean)

	if(min(c(f,m)) <= 1) {print('not enough individuals')} else {

	design = model.matrix.lm(~ 0 + sex + age + PMI + RIN + nreads + total_cells, data=meta_data_now_sub, na.action="na.pass")
	
	f = length(subset(meta_data_now_sub, sex == 'Female')$sex)
	m = length(subset(meta_data_now_sub, sex == 'Male')$sex)
	
	if(min(c(f,m)) <= 1) {print('not enough individuals')} else {
	
	y = voomWithQualityWeights(d0, plot = F)	
	
	if(dim(design)[1] <= dim(design)[2]) {print('not enough individuals')} else {
	
	fit = lmFit(y, design)
	contr = makeContrasts(sexMale - sexFemale, levels = colnames(coef(fit)))
	tmp = contrasts.fit(fit, contr)
	tmp = eBayes(tmp)
	res = data.frame(beta = tmp$coefficients, s2 = tmp$s2.post, stdev = tmp$stdev.unscaled, p = tmp$p.value)
	colnames(res) = c('beta','s2','stdev','p')	
	res$se = sqrt(res$s2) * res$stdev
	colnames(res) = c('log2FoldChange','s2','stdev','pvalue','lfcSE')
	
	saveRDS(res, file = paste(region_now,celltype_now,meta_data_now$new.id[i],"loo_SEX_DE_res.rds",sep="_"))

	TRADE_output = TRADE(mode = "univariate",
                     results1 = res,
                     annot_table = NULL,
                     genes_exclude = NULL,
                     n_sample = NULL)
                     
	out[i,1] = TRADE_output$distribution_summary$transcriptome_wide_impact

	}}}}
	
rownames(out) = meta_data_now$new.id
saveRDS(out, file = paste(region_now,celltype_now,"DE-loo-TRADE-jk.rds",sep="_"))

obs = readRDS(paste(region_now,celltype_now,'TRADE_sex.rds',sep = "_"))
obs = obs$distribution_summary$transcriptome_wide_impact
jkmean = mean(out$V1)
jkssq = sum((out$V1 - jkmean)^2)
n = length(meta_data_now$new.id)
jkse = sqrt((n - 1) / n * jkssq)
t_statistic = obs / jkse
p_value = 2 * pt(-abs(t_statistic), df = n - 1)

saveRDS(p_value, file = paste(region_now,celltype_now,"TRADE-jk-pval.rds",sep="_"))

print('done')
