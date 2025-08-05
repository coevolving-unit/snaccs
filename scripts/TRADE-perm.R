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
library(DESeq2)
library(TRADE)

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
meta_data_now = unique(meta_data_now[,c('new.id','region','RIN','PMI','individual',
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

# normalization

print('normalization')

d0 = DGEList(pseudo)
d0 = d0[keep, , keep.lib.sizes=FALSE]
d0 = calcNormFactors(d0, method = 'TMM')
dds = DESeqDataSetFromMatrix(countData = d0$counts,
                              colData = meta_data_now,
                              design = ~ sex)
dds = DESeq(dds)
cooks = assays(dds)[["cooks"]]
cooks_outlier = cooks > qf(0.99, df1 = 2, df2 = ncol(dds) - 2)
outlier_genes = rownames(dds)[rowSums(cooks_outlier) > 0]
dds_clean = replaceOutliers(dds)
d0$counts = counts(dds_clean)

y = voomWithQualityWeights(d0, plot = F)
	
f = length(subset(meta_data_now, sex == 'Female')$sex)
m = length(subset(meta_data_now, sex == 'Male')$sex)

print('scale all continuous variables')

meta_data_now$RIN = scale(meta_data_now$RIN)
meta_data_now$PMI = scale(meta_data_now$PMI)
meta_data_now$nreads = scale(meta_data_now$nreads)
meta_data_now$total_cells = scale(meta_data_now$total_cells)
meta_data_now$age = scale(meta_data_now$age)
meta_data_now$sex = as.factor(meta_data_now$sex)
meta_data_now$batch = as.factor(meta_data_now$batch)

if(min(c(f,m)) <= 1) {print('not enough individuals')} else {

	design = model.matrix.lm(~ 0 + sex + batch + age + PMI + RIN + nreads + total_cells, data=meta_data_now, na.action="na.pass")
	
	f = length(subset(meta_data_now, sex == 'Female')$sex)
	m = length(subset(meta_data_now, sex == 'Male')$sex)
	
	if(min(c(f,m)) <= 1) {print('not enough individuals')} else {
				
	if(dim(design)[1] <= dim(design)[2]) {print('not enough individuals')} else {
	
	comp = data.frame()
	out = data.frame()
	for(k in 1:1000){
		
		meta_perm = meta_data_now
		meta_perm$sex = sample(meta_perm$sex)
		
		design = model.matrix.lm(~ 0 + sex + batch + age + PMI + RIN + nreads + total_cells, data=meta_perm, na.action="na.pass")
		
		fit = lmFit(y, design)
		contr = makeContrasts(sexMale - sexFemale, levels = colnames(coef(fit)))
		tmp = contrasts.fit(fit, contr)
		tmp = eBayes(tmp)
		res = data.frame(beta = tmp$coefficients, s2 = tmp$s2.post, stdev = tmp$stdev.unscaled, p = tmp$p.value)
		colnames(res) = c('beta','s2','stdev','p')	
		res$se = sqrt(res$s2) * res$stdev
		colnames(res) = c('log2FoldChange','s2','stdev','pvalue','lfcSE')
	
		if(nrow(comp)==0) {comp = data.frame(tmp$coefficients)} else {
		
		comp[,k] = tmp$coefficients}
		
		TRADE_output = TRADE(mode = "univariate",
                     results1 = res,
                     annot_table = NULL,
                     genes_exclude = NULL,
                     n_sample = NULL)
                     
		out[k,1] = TRADE_output$distribution_summary$transcriptome_wide_impact
	
		}
	
	print('comparing to original DE')
	
	original = readRDS(paste(region_now,celltype_now,"SEX_DE_limma_outliers_corrected.rds",sep="_"))
	
	beta_perm = data.frame()
	for(k in 1:length(rownames(original))){
			obs = original$beta[k]
			co = comp[rownames(original)[k],]
			p_value = sum(abs(co) > abs(obs))/1000
			beta_perm[k,1] = rownames(original)[k]
			beta_perm[k,2] = p_value
			}
		
	colnames(beta_perm) = c('gene','p')
	beta_perm$padj = p.adjust(beta_perm$p)
	saveRDS(beta_perm, file = paste(region_now,celltype_now,"DE_perm.rds",sep="_"))
	
	print('comparing to original TRADE')
		
	obs = readRDS(paste(region_now,celltype_now,'TRADE_sex_outliers.rds',sep = "_"))
	obs = obs$distribution_summary$transcriptome_wide_impact
	out$p = sum(out$V1 > obs, na.rm = T)/1000
	
	saveRDS(out, file = paste(region_now,celltype_now,"TRADE_perm.rds",sep="_"))
	

}}}

print('done')
