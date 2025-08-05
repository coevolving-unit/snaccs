print('loading data')

library(mashr)
library(maditr)

`%!in%` = Negate(`%in%`)

all_res = readRDS('all_res_outliers_SEX_class_corrected.rds')

meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,41:151)])
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
all_res$key = paste(all_res$region, all_res$celltype)
all_res = subset(all_res, key %!in% remove$key)

all_res$beta = as.numeric(all_res$beta)
Bhat = maditr::dcast(data = all_res, formula = gene ~ celltype + region, value.var = 'beta', sep = "-")
Bhat = data.frame(Bhat)
rownames(Bhat) = Bhat$gene
Bhat = Bhat[,-1, drop=FALSE]

all_res$se = sqrt(all_res$s2) * all_res$stdev 

all_res$se = as.numeric(all_res$se)
Shat = maditr::dcast(data = all_res, formula = gene ~ celltype + region, value.var = 'se', sep = "-")
Shat = data.frame(Shat)
rownames(Shat) = Shat$gene
Shat = Shat[,-1, drop=FALSE]

Bhat = as.matrix(Bhat)
Shat = as.matrix(Shat)

mash.data = mash_set_data(Bhat,Shat)

set.seed(888)

U.c = cov_canonical(mash.data)

print('mash_1by1')

m.1by1 = mash_1by1(mash.data)
mash.thresh = 0.2
strong.subset = get_significant_results(m.1by1, thresh = mash.thresh)
data.strong = mash_set_data(Bhat[strong.subset,],Shat[strong.subset,])

## approach 2 - 19 PCs + EM null

print('data driven covariance')

U.pca = cov_pca(data.strong,19)
U.ed = cov_ed(data.strong, U.pca, maxiter = 10)
U.c = cov_canonical(mash.random)

print('mash_estimate_corr_em')

V.em = mash_estimate_corr_em(mash.data, c(U.c, U.ed), details = TRUE, max_iter = 5)
model = V.em$mash.model

idx <- which(is.na(Bhat))
model$result$PosteriorMean[idx] <- NA
model$result$PosteriorSD[idx] <- NA
model$result$NegativeProb[idx] <- NA
model$result$lfsr[idx] <- NA

print('save results')

saveRDS(model, file = 'mashr_sex_class_approach2.rds') 

print('done')
