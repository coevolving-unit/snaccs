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
library(sva)

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

	design = model.matrix(~ 0 + sex + age + PMI + RIN + nreads + total_cells, data=meta_data_now, na.action="na.pass")
	
	f = length(subset(meta_data_now, sex == 'Female')$sex)
	m = length(subset(meta_data_now, sex == 'Male')$sex)
	
	if(min(c(f,m)) <= 1) {print('not enough individuals')} else {
	
	write.csv(table(meta_data_now$sex), file = paste(region_now,celltype_now,"sex_N.csv",sep="_"))
		
	y = voomWithQualityWeights(d0, plot = F)
	
	print('saving normalized expression matrix')
	
	saveRDS(y$E, file = paste(region_now,celltype_now,"normalized_expression_corrected.rds",sep="_"))
	
	print('clustering samples')
	
	M = cor(y$E, use = 'pairwise.complete.obs', method = 'pearson')
	d = dist(1-M)
	cl = hclust(d)
	
	pdf(file = paste(region_now,celltype_now,"hclust_corrected.pdf",sep="_"))
	plot(cl)
	dev.off()

	# differential expression (sex)
	
	print('differential expression analysis')
	
	if(dim(design)[1] <= dim(design)[2]) {print('not enough individuals')} else {
	
	fit = lmFit(y, design)
	contr = makeContrasts(sexMale - sexFemale, levels = colnames(coef(fit)))
	tmp = contrasts.fit(fit, contr)
	tmp = eBayes(tmp)
	head(tmp$coefficients[order(tmp$coefficients[,1]),])
	tail(tmp$coefficients[order(tmp$coefficients[,1]),])

	print('saving DE output')

	saveRDS(tmp, file = paste(region_now,celltype_now,"fit_sex_corrected.rds",sep="_"))
	
	res = data.frame(beta = tmp$coefficients, s2 = tmp$s2.post, stdev = tmp$stdev.unscaled, p = tmp$p.value)
	colnames(res) = c('beta','s2','stdev','p')

	saveRDS(res, file = paste(region_now,celltype_now,"SEX_DE_res_corrected.rds",sep="_"))
	
	# extracting age results
	
	print('extracting age effects')

	fit2 = eBayes(fit)
	saveRDS(fit2, file = paste(region_now,celltype_now,"fit_age.rds",sep="_"))
	res = data.frame(beta = fit2$coefficients[,'age'], s2 = fit2$s2.post, 
					stdev = fit2$stdev.unscaled[,'age'], p = fit2$p.value[,'age'])
	colnames(res) = c('beta','s2','stdev','p')
	saveRDS(res, file = paste(region_now,celltype_now,"AGE_DE_res_corrected.rds",sep="_"))
	
	print('DESeq + replace outliers')
	
	dds = DESeqDataSetFromMatrix(countData = d0$counts,
                              colData = meta_data_now,
                              design = ~ sex)
	dds = DESeq(dds)
	cooks = assays(dds)[["cooks"]]
	cooks_outlier = cooks > qf(0.99, df1 = 2, df2 = ncol(dds) - 2)
	outlier_genes = rownames(dds)[rowSums(cooks_outlier) > 0]
	dds_clean = replaceOutliers(dds)
	design(dds_clean) = ~ 0 + sex + age + PMI + RIN + nreads + total_cells
	dds_clean = DESeq(dds_clean)
	res2 = results(dds_clean, contrast=c("sex","Male","Female"))
	saveRDS(res2, file = paste(region_now,celltype_now,"SEX_DE_deseq2_outliers_corrected.rds",sep="_"))

	print('voom limma without outliers')
	
	d0$counts = counts(dds_clean)
	y = voomWithQualityWeights(d0, plot = F)
	
	print('saving normalized expression matrix')
	saveRDS(y$E, file = paste(region_now,celltype_now,"normalized_expression_corrected_outliers.rds",sep="_"))

	fit = lmFit(y, design)
	contr = makeContrasts(sexMale - sexFemale, levels = colnames(coef(fit)))
	tmp = contrasts.fit(fit, contr)
	tmp = eBayes(tmp)
	res = data.frame(beta = tmp$coefficients, s2 = tmp$s2.post, stdev = tmp$stdev.unscaled, p = tmp$p.value)
	colnames(res) = c('beta','s2','stdev','p')
	saveRDS(res, file = paste(region_now,celltype_now,"SEX_DE_limma_outliers_corrected.rds",sep="_"))
	
	fit2 = eBayes(fit)
	saveRDS(fit2, file = paste(region_now,celltype_now,"fit_model_outliers.rds",sep="_"))
	res = data.frame(beta = fit2$coefficients[,'age'], s2 = fit2$s2.post, 
					stdev = fit2$stdev.unscaled[,'age'], p = fit2$p.value[,'age'])
	colnames(res) = c('beta','s2','stdev','p')
	saveRDS(res, file = paste(region_now,celltype_now,"AGE_DE_limma_outliers_corrected.rds",sep="_"))
	
	print('sva and limma')
		
	mod = model.matrix(~ 0 + sex + age + PMI + RIN + nreads + total_cells, data=meta_data_now, na.action="na.pass")
	mod0 = model.matrix(~ 0 + age + PMI + RIN + nreads + total_cells, data=meta_data_now, na.action="na.pass")
	
	n.sv = num.sv(d0$counts,mod,method="leek")
	n.sv = n.sv -1
	svobj = sva(d0$counts,mod,mod0,n.sv=n.sv)
	modSv = cbind(mod,svobj$sv)
	fit = lmFit(d0$counts,modSv)
	contr = makeContrasts(sexMale - sexFemale, levels = c(colnames(coef(fit))[which(colnames(coef(fit)) != "")],paste("SV",c(1:n.sv),sep="")))
	tmp = contrasts.fit(fit, contr)
	tmp = eBayes(tmp)
	res = data.frame(beta = tmp$coefficients, 
						s2 = tmp$s2.post, 
						stdev = tmp$stdev.unscaled, 
						p = tmp$p.value,
						nsv = n.sv)
	colnames(res) = c('beta','s2','stdev','p','nsv')

	saveRDS(res, file = paste(region_now,celltype_now,"SEX_DE_sva_limma_corrected.rds",sep="_"))

	print('DESeq + replace outliers + sva')
	
	ddssva = dds_clean
	colData(ddssva) = cbind(colData(ddssva), svobj$sv)
	v = paste("V",c(1:dim(svobj$sv)[2]),"+",sep="")
	vall = paste(v, collapse="")
	vall = str_sub(vall, 1, -2)
	des = paste("~ 0 + sex + age + PMI + RIN + nreads + total_cells +", vall, sep="")
	des_formula <- as.formula(des)
	design(ddssva) <- des_formula
	ddssva = DESeq(ddssva)
	ressva = results(ddssva, contrast=c("sex","Male","Female"))
	saveRDS(ressva, file = paste(region_now,celltype_now,"SEX_DE_deseq2_outliers_sva_corrected.rds",sep="_"))
	

}}}

print('done')
