library(tidyr)
library(data.table)
library(mashr)
library(pheatmap)
library(stringr)
library(ggplot2)
library(egg)
library(RColorBrewer)
library(dplyr)
library(NbClust)
library(ggpubr)
library(ggrepel)

#### compare different mashr results ####

sex1 = readRDS('/mashr/mashr_sex_class_approach1.rds')
sexbeta1 = get_pm(sex1)
sexlfsr1 = get_lfsr(sex1)
sexbeta1 = reshape1::melt(sexbeta1)
sexlfsr1 = reshape1::melt(sexlfsr1)
sexres1 = cbind(sexbeta1, sexlfsr1)
d = str_split_fixed(sexres1$Var1, pattern = '\\.', 2)
sexres1 = cbind(sexres1, d)
sexres1 = sexres1[,c(1,3,6:8)]
colnames(sexres1) = c('gene','beta','lfsr','cell','region')

sex2 = readRDS('/mashr/mashr_sex_class_approach2.rds')
sexbeta2 = get_pm(sex2)
sexlfsr2 = get_lfsr(sex2)
sexbeta2 = reshape2::melt(sexbeta2)
sexlfsr2 = reshape2::melt(sexlfsr2)
sexres2 = cbind(sexbeta2, sexlfsr2)
d = str_split_fixed(sexres2$Var2, pattern = '\\.', 2)
sexres2 = cbind(sexres2, d)
sexres2 = sexres2[,c(1,3,6:8)]
colnames(sexres2) = c('gene','beta','lfsr','cell','region')

sex3 = readRDS('/mashr/mashr_sex_class_approach3.rds')
sexbeta3 = get_pm(sex3)
sexlfsr3 = get_lfsr(sex3)
sexbeta3 = reshape2::melt(sexbeta3)
sexlfsr3 = reshape2::melt(sexlfsr3)
sexres3 = cbind(sexbeta3, sexlfsr3)
d = str_split_fixed(sexres3$Var2, pattern = '\\.', 2)
sexres3 = cbind(sexres3, d)
sexres3 = sexres3[,c(1,3,6:8)]
colnames(sexres3) = c('gene','beta','lfsr','cell','region')

sex4 = readRDS('/mashr/mashr_sex_class_approach4.rds')
sexbeta4 = get_pm(sex4)
sexlfsr4 = get_lfsr(sex4)
sexbeta4 = reshape2::melt(sexbeta4)
sexlfsr4 = reshape2::melt(sexlfsr4)
sexres4 = cbind(sexbeta4, sexlfsr4)
d = str_split_fixed(sexres4$Var2, pattern = '\\.', 2)
sexres4 = cbind(sexres4, d)
sexres4 = sexres4[,c(1,3,6:8)]
colnames(sexres4) = c('gene','beta','lfsr','cell','region')

sex5 = readRDS('/mashr/mashr_sex_class_approach5.rds')
sexbeta5 = get_pm(sex5)
sexlfsr5 = get_lfsr(sex5)
sexbeta5 = reshape2::melt(sexbeta5)
sexlfsr5 = reshape2::melt(sexlfsr5)
sexres5 = cbind(sexbeta5, sexlfsr5)
d = str_split_fixed(sexres5$Var2, pattern = '\\.', 2)
sexres5 = cbind(sexres5, d)
sexres5 = sexres5[,c(1,3,6:8)]
colnames(sexres5) = c('gene','beta','lfsr','cell','region')

# compare likelihoods

get_loglik(sex1) # -11590304 # 4
get_loglik(sex2) # -11299407 # 2
get_loglik(sex3) # -11299200 # 1
get_loglik(sex4) # -11646309 # 5
get_loglik(sex5) # -11301376 # 3

# estimated mixture proportions

View(data.frame(round(get_estimated_pi(sex1),3))) # ED_tPCA + ED_PCA_7 + equal + ED_PCA_6 + OPC ILTO + ED_PCA_7 + Micro FG
View(data.frame(round(get_estimated_pi(sex2),3))) # OPC ILTO + equal + ED_PCA_7 + ED_tPCA + Micro FG
View(data.frame(round(get_estimated_pi(sex3),3))) # OPC ILTO + equal + ED_tPCA + ED_PCA_7 + Micro FG
View(data.frame(round(get_estimated_pi(sex4),3))) # mostly equal effects + OPC ILTO + het2 + Micro FG
View(data.frame(round(get_estimated_pi(sex5),3))) # OPC ILTO + equal + het2 + het2 + Micro FG

pl = data.frame(mi = get_estimated_pi(sex3))
pl$group = rownames(pl)

ggplot(pl, aes(x = reorder(group, -mi), y = mi)) +
  geom_bar(stat = 'identity') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

# compute proportion of shared signals across conditions

get_pairwise_sharing(sex1)
get_pairwise_sharing(sex2)
get_pairwise_sharing(sex3)
get_pairwise_sharing(sex4)
get_pairwise_sharing(sex5)

# view individual covariances

corrplot::corrplot(sex3$fitted_g$Ulist$ED_PCA_7, 
                   order = 'hclust',
                   cl.cex = 0.2, tl.cex = 0.4,
                   tl.col = 'black')

corrplot::corrplot(sex3$fitted_g$Ulist$ED_tPCA, 
                   order = 'hclust',
                   cl.cex = 0.2, tl.cex = 0.4,
                   tl.col = 'black')

corrplot::corrplot(sex3$fitted_g$Ulist$ED_PCA_6, 
                   order = 'hclust',
                   cl.cex = 0.2, tl.cex = 0.4,
                   tl.col = 'black')

corrplot::corrplot(sex3$fitted_g$Ulist$ED_FLASH_non_neg_4, 
                   order = 'hclust',
                   cl.cex = 0.2, tl.cex = 0.4,
                   tl.col = 'black')

corrplot::corrplot(sex3$fitted_g$Ulist$ED_tFLASH_non_neg, 
                   order = 'hclust',
                   cl.cex = 0.2, tl.cex = 0.4,
                   tl.col = 'black')

corrplot::corrplot(sex3$fitted_g$Ulist$ED_tFLASH_non_neg, 
                   order = 'hclust',
                   cl.cex = 0.2, tl.cex = 0.4,
                   tl.col = 'black')

#### end ####

#### compare pre- and post-mashr ####

# mashr

mash = readRDS('/mashr/mashr_sex_class_approach3.rds')

sexbeta = get_pm(mash)
sexlfsr = get_lfsr(mash)
sexbeta2 = reshape2::melt(sexbeta)
sexlfsr2 = reshape2::melt(sexlfsr)
sexres = cbind(sexbeta2, sexlfsr2)
d = str_split_fixed(sexres$Var2, pattern = '\\.', 2)
sexres = cbind(sexres, d)
sexres = sexres[,c(1,3,6:8)]
colnames(sexres) = c('gene','beta','lfsr','cell','region')
#sexres$bias = ifelse(sexres$beta > 0 & sexres$lfsr < 0.05, 'male', 'ns')
#sexres$bias = ifelse(sexres$beta < 0 & sexres$lfsr < 0.05, 'female', sexres$bias)
sexres$key = paste(sexres$gene, sexres$region, sexres$cell)

# pre-mashr

limma = readRDS('/mashr/all_res_outliers_SEX_class_corrected.rds')
colnames(limma)[8] = 'lfsr'
colnames(limma)[7] = 'cell'
limma$bias = ifelse(limma$beta > 0 & limma$lfsr < 0.05, 'male', 'ns')
limma$bias = ifelse(limma$beta < 0 & limma$lfsr < 0.05, 'female', limma$bias)
table(limma$bias, limma$region, limma$cell)
limma$key = paste(limma$gene, limma$region, limma$cell)

limma$se = sqrt(limma$s2) * limma$stdev 

ggplot(limma, aes(x = beta)) +
  geom_histogram() +
  facet_grid(cell ~ region, scales = 'free_y') +
  xlim(c(-2,2)) +
  theme_article()
ggplot(limma, aes(x = p)) +
  geom_histogram() +
  xlim(c(-0.1,1)) +
  facet_grid(cell ~ region, scales = 'free_y') +
  theme_article()
ggplot(limma, aes(x = se)) +
  geom_histogram() +
  xlim(c(-0.1,2)) +
  facet_grid(cell ~ region, scales = 'free_y') +
  theme_article()
ggplot(limma, aes(x = se, y = p)) +
  geom_point() +
  facet_grid(cell ~ region) +
  theme_article()
ggplot(limma, aes(x = se, y = beta)) +
  geom_point() +
  ylim(c(-1,1)) +
  facet_grid(cell ~ region) +
  theme_article()

# combine 

combo = merge(sexres, limma, by = 'key')
table(combo$region.x, combo$cell.x)
table(sign(combo$beta.x) == sign(combo$beta.y))
combo$match = ifelse(sign(combo$beta.x) == sign(combo$beta.y), 'match', 'diff')
diffsign = subset(combo, match == 'diff')
diffsign$key = paste(diffsign$region.x, diffsign$cell.x, diffsign$gene.x)

mashsig = subset(combo, lfsr.x < 0.05)
table(sign(mashsig$beta.x) == sign(mashsig$beta.y))
6526/(52605+6526) # 11% significant signals have opposite sign in limma

# plot # DEGs

mp = combo %>% group_by(region.x, cell.x, bias.x) %>% summarise(n = n())
lp = combo %>% group_by(region.x, cell.x, bias.y) %>% summarise(n = n())
cp = cbind(mp, lp$n)
colnames(cp) = c('region','celltype','bias','mashr','voom.limma')
cp = subset(cp, bias != 'ns')
cp$celltype = factor(cp$celltype, levels= c(exn.levels, inn.levels, other.levels))
cp$celltype = droplevels(cp$celltype)

ggplot(cp, aes(x = mashr, y = voom.limma)) + 
  geom_point(size = NA) +
  geom_smooth(method = 'lm') +
  geom_point(aes(color = celltype, shape = bias)) + 
  xlab('MASH') + ylab('voom-limma') +
  stat_cor(method = 'spearman', label.y.npc = 'top') +
  scale_color_manual(values = c(exn.colors, inn.colors, other.colors)) +
  theme_article()

cor.test(cp$mashr, cp$voom.limma, method = 'spearman')

ggplot(cp, aes(x = mashr, y = log(voom.limma))) + 
  geom_point(size = NA) +
  geom_smooth(method = 'lm') +
  geom_point(aes(color = celltype, shape = bias)) + 
  stat_cor(method = 'spearman', label.y.npc = 'top') +
  xlab('MASH') + ylab('voom-limma') +
  scale_color_manual(values = c(exn.colors, inn.colors, other.colors)) +
  theme_article()

cor.test(cp$mashr, log(cp$voom.limma), method = 'spearman')

cp_filtered <- subset(cp, !(region == "Inferior_Lateral_TemporalOccipital" & celltype == "OPC"))

ggplot(cp_filtered, aes(x = mashr, y = voom.limma)) + 
  geom_point(size = NA) +
  geom_smooth(method = 'lm') +
  geom_point(aes(color = celltype, shape = bias)) + 
  stat_cor(method = 'spearman', label.y.npc = 'top') +
  xlab('MASH') + ylab('voom-limma') +
  scale_color_manual(values = c(exn.colors, inn.colors, other.colors)) +
  theme_article()

cor.test(cp_filtered$mashr, cp_filtered$voom.limma, method = 'spearman')

# plot betas

combo$cell.x = factor(combo$cell.x, levels= c(exn.levels, inn.levels, other.levels))
combo$cell.x = droplevels(combo$cell.x)

ggplot(combo, aes(x = beta.x, y = beta.y)) + 
  geom_point(size = NA) +
  geom_smooth(method = 'lm') +
  geom_point() + 
  facet_grid(cell.x~region.x) +
  stat_cor(method = 'spearman') +
  xlab('MASH') + ylab('voom-limma') +
  theme_article()

cor.test(combo$beta.x, combo$beta.y, method = 'spearman')

cf = subset(combo, gene.x %!in% ygenes$gene)
cf = subset(cf, gene.x != 'XIST')

ggplot(cf, aes(x = beta.x, y = beta.y)) + 
  geom_point(size = NA) +
  geom_smooth(method = 'lm') +
  geom_point() + 
  facet_grid(cell.x~region.x) +
  stat_cor(method = 'spearman', label.y.npc = 'top') +
  xlab('MASH') + ylab('voom-limma') + 
  theme_article()

cor.test(cf$beta.x, cf$beta.y, method = 'spearman')

# find flips

sig_mash = subset(combo, lfsr.x < 0.05)
length(unique(sig_mash$gene.x))
dim(sig_mash)
table(sig_mash$beta.x > 0, sig_mash$beta.y > 0) 
# 6526 flip (from ns in limma to sig in mash) / 59131 total sig mash = 11%

sig_both = subset(combo, lfsr.x < 0.05 & lfsr.y < 0.05)
dim(sig_both)
table(sig_both$beta.x > 0, sig_both$beta.y > 0) 
# no genes are sig limma and then flip signs to be sig mashr

#### end ####

#### plot mashr sex effects (volcano + bar plots) ####

sex = readRDS('/Users/decasienar/Desktop/Papers/snRNAseq/alevin_mapped_20k/mashr_corrected/mashr_sex_class_approach3.rds')

sexbeta = get_pm(sex)
sexlfsr = get_lfsr(sex)
dim(sexbeta)

sexbeta2 = reshape2::melt(sexbeta)
sexlfsr2 = reshape2::melt(sexlfsr)
sexres = cbind(sexbeta2, sexlfsr2)
d = str_split_fixed(sexres$Var2, pattern = '\\.', 2)
sexres = cbind(sexres, d)
sexres = sexres[,c(1,3,6:8)]
colnames(sexres) = c('gene','beta','lfsr','cell','region')
sexres$key = paste(sexres$region, sexres$cell, sexres$gene)
#sexres$beta = ifelse(sexres$key %in% diffsign$key, 0, sexres$beta) # keep mashr as is
#sexres$lfsr = ifelse(sexres$key %in% diffsign$key, 1, sexres$lfsr) # keep mashr as is
sexres$bias = ifelse(sexres$beta > 0 & sexres$lfsr < 0.05, 'male', 'ns')
sexres$bias = ifelse(sexres$beta < 0 & sexres$lfsr < 0.05, 'female', sexres$bias)
table(sexres$bias)
sexres$key = paste(sexres$region, sexres$cell)

meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,44:154)]) 
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var2, check$Var2)
remove = subset(check, Freq < 10)
sexres = subset(sexres, key %!in% remove$key)

saveRDS(sexres, file = 'sexres.rds')

sexres = readRDS('sexres.rds')
bb = sexres %>% group_by(gene, bias) %>% summarise(n=n())
colnames(bb)[1] = 'external_gene_name'
gene2chrom = readRDS('gene2chrom_mash.rds')
bb = merge(bb, gene2chrom, by = 'external_gene_name', all.x = T)
View(subset(bb, n == 114 & bias != 'ns'))
length(unique((subset(bb, n == 114 & bias != 'ns')$external_gene_name)))
length(unique((subset(bb, n == 114 & bias != 'ns' & chromosome_name %!in% c('X','Y'))$external_gene_name)))
length(unique((subset(bb, n == 114 & bias == 'male' & chromosome_name %!in% c('X','Y'))$external_gene_name)))
length(unique((subset(bb, n == 114 & bias == 'female' & chromosome_name %!in% c('X','Y'))$external_gene_name)))
c = as.character(subset(bb, n == 114 & bias != 'ns')$external_gene_name)
table(c)

length(unique(subset(sexres, bias %in% c('male','female'))$gene))
length(unique(subset(bb, bias %in% c('male','female') & chromosome_name %in% c(1:22))$external_gene_name))

pl = sexres %>% group_by(gene, bias) %>% summarise(n = n())
pl = subset(pl, bias %in% c('male','female'))
length(subset(pl, n == 114)$gene)

pl = sexres %>% group_by(gene, bias, region) %>% summarise(n = n())
pl = subset(pl, bias %in% c('male','female'))
pl = unique(pl[,c('gene', 'region')])
pl = pl %>% group_by(gene) %>% summarise(n = n())
length(subset(pl, n == 6)$gene)
length(subset(pl, n == 1)$gene)
length(subset(pl, n > 1)$gene)

pl = sexres %>% group_by(gene, bias, cell) %>% summarise(n = n())
pl = subset(pl, bias %in% c('male','female'))
pl = unique(pl[,c('gene', 'cell')])
pl = pl %>% group_by(gene) %>% summarise(n = n())
length(subset(pl, n == 19)$gene)
length(subset(pl, n == 1)$gene)
length(subset(pl, n > 1)$gene)

View(sexres %>% group_by(region) %>% summarise(mean = mean(beta, na.rm = T), absmean = mean(abs(beta),na.rm = T)))
View(sexres %>% group_by(cell) %>% summarise(mean = mean(beta, na.rm = T), absmean = mean(abs(beta),na.rm = T)))
View(sexres %>% group_by(region, cell) %>% summarise(mean = mean(beta, na.rm = T), absmean = mean(abs(beta),na.rm = T)))

mod = aov(abs(beta) ~ region + cell, data = sexres)
options(digits = 16)
summary(mod)
summary(mod)[[1]][["Pr(>F)"]]
c=TukeyHSD(mod)
df = data.frame(rbind(c$region, c$cell))
write.csv(df, file = 'tukey-abs-beta-region-cell.csv')

b = sexres %>% group_by(gene, bias) %>% summarise(n = n())
b$key = paste(b$gene, b$bias)
max(b$n)

sexres_min <- sexres %>%
  group_by(gene) %>%
  slice_min(order_by = lfsr, with_ties = FALSE, na_rm = T) %>%
  ungroup()
dim(sexres_min)
sexres_min$key = paste(sexres_min$gene, sexres_min$bias)
sexres_min = merge(sexres_min, b, by = 'key', all.x = T)
table(sexres_min$lfsr < 0.05)
sexres_min$beta = ifelse(abs(sexres_min$beta) >=5, sign(sexres_min$beta)*5, sexres_min$beta)
sexres_min$lfsr = ifelse(sexres_min$lfsr <= 1e-20, 1e-20, sexres_min$lfsr)
sexres_min$n = ifelse(sexres_min$bias.x == 'ns', 1, sexres_min$n)
colnames(sexres_min)[2] = 'gene'

# hsap = useEnsembl(biomart = 'ensembl', dataset='hsapiens_gene_ensembl', version = 98) 
# gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), 
#                     filters = 'external_gene_name', values = unique(sexres$gene), mart = hsap)
# table(unique(sexres$gene) %in% gene2chrom$external_gene_name)
# saveRDS(gene2chrom, file = 'gene2chrom_mash.rds')

# hsap = useEnsembl(biomart = 'ensembl', dataset='hsapiens_gene_ensembl', version = 222)
# parcheck = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), 
#                                        filters = 'external_gene_name', values = unique(sexres$gene), mart = hsap)
# par = parcheck %>% group_by(external_gene_name) %>% summarise(n = n())
# par = subset(par, n == 2)
# par = subset(parcheck, external_gene_name %in% par$external_gene_name)
# par = subset(par, chromosome_name %in% c('X','Y'))
# par = unique(par[,c(2:3)])
# par = par %>% group_by(external_gene_name) %>% summarise(n = n())
# par = subset(par, n == 2)
# saveRDS(par, file = 'par.rds')

par = readRDS('par.rds')
par = rbind(par, data.frame(external_gene_name = 'PPP2R3B', n = 2))
gene2chrom = readRDS('gene2chrom_mash.rds')
table(gene2chrom$chromosome_name)
gene2chrom = subset(gene2chrom, chromosome_name %in% c(1:22, 'X', 'Y'))
table(table(gene2chrom$external_gene_name))
table(unique(sexres$gene) %in% gene2chrom$external_gene_name)
unique(sexres$gene)[which(unique(sexres$gene) %!in% gene2chrom$external_gene_name)]
gene2chrom[which(gene2chrom$external_gene_name %in% par$external_gene_name),'chromosome_name'] = 'PAR'
colnames(gene2chrom)[2] = 'gene'

sexres_min = merge(sexres_min, unique(gene2chrom[,c('gene','chromosome_name')]), by = 'gene', all.x = T)
table(sexres_min$chromosome_name)
sexres_min$chrom2 = ifelse(sexres_min$chromosome_name == 'X' | sexres_min$chromosome_name == 'Y' | sexres_min$chromosome_name == 'PAR', sexres_min$chromosome_name, 'Autosomal')
sexres_min = sexres_min[,-11]
sexres_min = unique(sexres_min)

jitter <- position_jitter(width = 0.3, height = 0.3)
sexres_min = sexres_min[complete.cases(sexres_min),]
ggplot(sexres_min) +
  geom_point(aes(y = -log10(lfsr), x = beta, size = n, shape = chrom2, 
                 fill = bias.x, color = bias.x), alpha = 0.4, position = jitter) +
  scale_fill_manual(values = c(sex.colors,"#666666")) +
  scale_color_manual(values = c(sex.colors,"#666666")) +
  scale_shape_manual(values = c(21, 22, 4, 24)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  #geom_text_repel(aes(y = -log10(lfsr), x = beta, label = gene.x), size = 3, max.overlaps = 30) +
  theme_article() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_blank(),
        legend.position = 'none')

genes_to_label <- c("CHM", "DDX3X", "JPX", "KDM5C", "KDM6A", "UBA1", "ZFX", 
                      "USP9X", "XIST", "LINC00894")
genes_to_label <- par$external_gene_name

ggplot(sexres_min) +
  geom_point(aes( y = -log10(lfsr), x = beta,size = n,shape = chrom2, fill = chrom2, color = chrom2), alpha = 0.4, position = position_jitter(width = 0.1, height = 0)) +
  #scale_fill_manual(values = c(sex.colors, "#666666")) +
  #scale_color_manual(values = c(sex.colors, "#666666")) +
  scale_shape_manual(values = c(21, 22, 4, 24)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_text_repel(
    data = subset(sexres_min, gene %in% genes_to_label),
    aes(y = -log10(lfsr), x = beta, label = gene),
    size = 3,
    max.overlaps = 30) +
  theme_article() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_blank(),
    legend.position = 'none')

pl = sexres[complete.cases(sexres),]
pl = merge(pl, unique(gene2chrom[,c('gene','chromosome_name')]), by = 'gene', all.x = T)
unique(pl[which(is.na(pl$chromosome_name)),'gene'])
pl = pl[complete.cases(pl),]
pl$chrom2 = ifelse(pl$chromosome_name == 'X' | pl$chromosome_name == 'Y' | pl$chromosome_name == 'PAR', pl$chromosome_name, 'Autosomal')
pl = pl %>% group_by(cell, region, bias, chrom2) %>% summarise(n = n())
pl$key = paste(pl$cell, pl$region)
pl2 = pl %>% group_by(cell, region) %>% summarise(sum = sum(n))
pl2$key = paste(pl2$cell, pl2$region)
pl = merge(pl, pl2, by = 'key')
pl$prop = pl$n / pl$sum
pl$cell.x = factor(pl$cell.x,levels = c(exn.levels, inn.levels, other.levels))
pl$region.x = factor(pl$region.x)
levels(pl$region.x) = regions_short

ggplot(pl) +
  geom_bar(stat = 'identity', aes(x = cell.x, y = prop, fill = bias)) +
  scale_fill_manual(values = c(sex.colors,'#666666')) +
  facet_wrap(~region.x, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

pl = subset(pl, bias != 'ns')

ggplot(pl) +
  geom_bar(stat = 'identity', aes(x = cell.x, y = n, fill = bias)) +
  scale_fill_manual(values = c(sex.colors)) +
  facet_wrap(~region.x, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

ggplot(pl) +
  geom_bar(stat = 'identity', aes(x = cell.x, y = n, fill = bias)) +
  scale_fill_manual(values = c(sex.colors)) +
  facet_grid(chrom2~region.x, scales = 'free') +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 10),
        axis.title = element_blank())

View(pl %>% group_by(cell.x) %>% summarise(sum = sum(n), mean = mean(n)))
View(pl %>% group_by(region.x) %>% summarise(sum = sum(n), mean = mean(n)))
View(pl %>% group_by(cell.x, region.x) %>% summarise(sum = sum(n), mean = mean(n)))

degs = pl %>% group_by(cell.x, region.x) %>% summarise(prop = sum(prop), n = sum(n))
colnames(degs) = c('celltype','region','prop','n')
degs$rankdeg = rank(degs$n, ties.method = "random")  
degs$rankprop = rank(degs$prop, ties.method = "random")  

mod = aov(prop ~ celltype + region, degs)
summary(mod)
c=TukeyHSD(mod)
df = data.frame(rbind(c$region, c$cell))
write.csv(df, file = 'tukey-prop-sexbiased-region-cell.csv')

pl = subset(sexres, bias %in% c('male','female'))
pl$cell = factor(pl$cell,levels = c(exn.levels, inn.levels, other.levels))
pl$region = factor(pl$region)
levels(pl$region) = regions_short

ygenes = subset(gene2chrom, chromosome_name == 'Y')
xgenes = subset(gene2chrom, chromosome_name == 'X')

pl = subset(pl, gene %!in% c('XIST','TSIX'))
pl = subset(pl, gene %!in% ygenes$external_gene_name) # x and auto
pl = subset(pl, gene %in% xgenes$external_gene_name) # x only
pl = subset(pl, gene %!in% c(ygenes$external_gene_name,xgenes$external_gene_name)) # auto only

ggplot(pl, aes(x = cell, y = beta)) +
  geom_violin(data = subset(pl, beta > 0), aes(x = cell, y = beta), fill = sex.colors[2], width = 2) +
  geom_violin(data = subset(pl, beta < 0), aes(x = cell, y = beta), fill = sex.colors[1], width = 2) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

ggplot(pl, aes(x = abs(beta))) +
  geom_density() +
  facet_grid(region~cell) +
  theme_article() +
  scale_x_continuous(limits = c(0,2), breaks = c(0,1)) +
  theme(legend.position = 'none',
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_blank())

#### end ####

#### plot mashr sex effects (correlation/block modelling) #### 

sexbetanow = sexbeta

# autosomal only
sexbetanow = sexbetanow[which(rownames(sexbetanow) %in% subset(gene2chrom, chromosome_name %in% c(1:22))$gene),]

# correlation
sexc = WGCNA::cor(sexbetanow, method = 'spearman', use = 'pairwise.complete.obs')
saveRDS(sexc, file = 'sexc.rds')

# format
rownames(sexc)=gsub("Angular_Gyrus",'AngGyr', rownames(sexc))
rownames(sexc)=gsub("Caudal_Insula",'CauIns', rownames(sexc))
rownames(sexc)=gsub("Fusiform_Gyrus",'FusGyr', rownames(sexc))
rownames(sexc)=gsub("Intraparietal_sulcus",'IntParS', rownames(sexc))
rownames(sexc)=gsub("Inferior_Lateral_TemporalOccipital",'InfLTO', rownames(sexc))
rownames(sexc)=gsub("Retrosplenial_cortex",'RetCor', rownames(sexc))
colnames(sexc)=gsub("Angular_Gyrus",'AngGyr', colnames(sexc))
colnames(sexc)=gsub("Caudal_Insula",'CauIns', colnames(sexc))
colnames(sexc)=gsub("Fusiform_Gyrus",'FusGyr', colnames(sexc))
colnames(sexc)=gsub("Intraparietal_sulcus",'IntParS', colnames(sexc))
colnames(sexc)=gsub("Inferior_Lateral_TemporalOccipital",'InfLTO', colnames(sexc))
colnames(sexc)=sub("Retrosplenial_cortex",'RetCor', colnames(sexc))

# block modelling

library(blockmodels)
set.seed('123') 
diff_matrix_op <- BM_gaussian("SBM_sym", sexc, explore_max = 20) 
diff_matrix_op$estimate()
diff_matrix_op = saveRDS(diff_matrix_op, file = 'diff_matrix_op.rds')

optimum_clusters <- which.max(diff_matrix_op$ICL) 
print(optimum_clusters)
diff_cluster_memberships <- apply(diff_matrix_op$memberships[[optimum_clusters]]$Z, 1, function(x) which.max(x)) 
table(diff_cluster_memberships)
o = data.frame(scale=colnames(sexc), sbm_cluster=diff_cluster_memberships) %>% arrange(sbm_cluster) #view as dataframe
#saveRDS(o, 'blocks.rds')

o = readRDS('blocks.rds')

diff_block.reorder <- c(1:optimum_clusters)[order(colMeans(diff_matrix_op$model_parameters[[optimum_clusters]]$mu))] # mu takes mean of upper triangle of each block
diff_cluster_reordered <- as.factor(diff_cluster_memberships)
levels(diff_cluster_reordered) <- diff_block.reorder # reorder the factor variable based on the correlation
diff_cluster_labels <- data.frame(scale = colnames(sexc), 
                                  sbm_cluster = diff_cluster_memberships, 
                                  sbm_cluster_ordered = factor(diff_cluster_memberships,levels = order(colMeans(diff_matrix_op$model_parameters[[optimum_clusters]]$mu)))) %>%
arrange(diff_cluster_reordered) # 3 7 2 6 5 4 1
diff_cluster_labels$names <- factor(diff_cluster_labels$sbm_cluster_ordered)
diff_cluster_reordered <- as.numeric(as.character(diff_cluster_reordered))
diff_cluster_order <- as.numeric(as.character(diff_cluster_labels$sbm_cluster))
custom_cols <- diff_cluster_labels$scale
custom_rows <- diff_cluster_labels$scale
col_indices <- match(custom_cols, colnames(sexc))
row_indices <- match(custom_rows, rownames(sexc))
matrix_ordered <- sexc[row_indices,col_indices]

library(igraph)
library(plotrix)

sbm_g_op_adj <- diff_matrix_op$model_parameters[[optimum_clusters]]$mu
adj_ordered <- sbm_g_op_adj[diff_block.reorder, diff_block.reorder]
View(adj_ordered)
adj_ordered[which(adj_ordered < 0.3)] = 0
diag(adj_ordered) = 0
mean(adj_ordered)
hist(adj_ordered)
colnames(adj_ordered) <- levels(diff_cluster_labels$names)
rownames(adj_ordered) <- levels(diff_cluster_labels$names)
row_max_values <- apply(adj_ordered, 1, max)
colnames_for_max <- apply(adj_ordered, 1, function(x) colnames(adj_ordered)[which.max(x)])
d = data.frame(from = rownames(adj_ordered),
               to = colnames_for_max,
               max = row_max_values)
View(d)
ordered_grph <- graph_from_adjacency_matrix(adj_ordered, 
                                            mode="upper", 
                                            weighted=T, 
                                            diag=F, add.colnames=T)
V(ordered_grph)$label <- colnames(adj_ordered)
V(ordered_grph)$label.dist <- 1
V(ordered_grph)$label.cex <- plotrix::rescale(colMeans(adj_ordered), c(1,1.5))
V(ordered_grph)$label.color <- "black"
V(ordered_grph)$size <- plotrix::rescale(colMeans(adj_ordered), c(5, 10))
E(ordered_grph)$width <- plotrix::rescale(get.edge.attribute(ordered_grph)$weight, c(0.5,10))
plot(ordered_grph)

cluster_df = data.frame(cluster = o$sbm_cluster)
cluster_df$cluster = as.factor(cluster_df$cluster)
rownames(cluster_df) = o$scale

ann_colors <- list(
  cluster = c("1" = "#F3C300", 
              "2" = "#875692",
              "3" = "#F38400",
              "4" = "#A1CAF1",
              "5" = "#BE0032",
              "6" = "#C2B280",
              "7" = "#848482" ,
              "8" = "#008856",
              "9" = "#E68FAC",
              "10" = "#0067A5",
              "11" = "darkblue",
              "12" = "darkgreen",
              "13" = "purple4")
)

pheatmap(matrix_ordered,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="PRGn")))(100),
         fontsize = 5,
         breaks = seq(-1,1, length.out=(101)), 
         clustering_method = 'average',
         cluster_rows=FALSE,
         border_color = NA,
         cluster_cols=FALSE,
         annotation_col = cluster_df,
         annotation_row = cluster_df,
         annotation_colors = ann_colors)

region_df = data.frame(sub(".*\\.", "", colnames(matrix_ordered)))
cell_df = data.frame(sub("\\..*", "", colnames(matrix_ordered)))
rownames(region_df) = rownames(cell_df) = rownames(matrix_ordered)
colnames(region_df) = 'region'
colnames(cell_df) = 'cell'

names(exn.colors) <- exn.levels
names(inn.colors) <- inn.levels
names(other.colors) <- other.levels

ncol = which(c(exn.levels, inn.levels, other.levels) %!in% cell_df$cell)

ann_colors <- list(
  region = c("AngGyr" = region.colors[1], 
             "CauIns" = region.colors[2],
             "FusGyr" = region.colors[3],
             "InfLTO" = region.colors[4],
             "IntParS" = region.colors[5],
             "RetCor" = region.colors[6]),
  cell = c(exn.colors, inn.colors, other.colors)[-c(ncol)]
)

col2 = c(class.colors[1],class.colors[1],class.colors[1],class.colors[1],class.colors[1],class.colors[1],class.colors[1],
         class.colors[2],class.colors[2],class.colors[2],class.colors[2],class.colors[2],class.colors[2],class.colors[2],
         class.colors[3],class.colors[3],class.colors[3],class.colors[3],class.colors[3])
names(col2)= names(c(exn.colors, inn.colors, other.colors)[-c(ncol)])

ann_colors <- list(
  region = c("AngGyr" = region.colors[1], 
             "CauIns" = region.colors[2],
             "FusGyr" = region.colors[3],
             "InfLTO" = region.colors[4],
             "IntParS" = region.colors[5],
             "RetCor" = region.colors[6]),
  cell = col2)

pheatmap(matrix_ordered,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="PRGn")))(100),
         fontsize = 5,
         breaks = seq(-1,1, length.out=(101)), 
         clustering_method = 'average',
         cluster_rows=FALSE,
         border_color = NA,
         cluster_cols=FALSE,
         annotation_col = region_df,
         annotation_row = cell_df,
         annotation_colors = ann_colors)

clus = unique(o$sbm_cluster)
cm = data.frame(row.names = rownames(sexbetanow))
for(i in 1:length(clus)){
  bn = sexbetanow[,which(colnames(sexbetanow) %in% subset(o, sbm_cluster == clus[i])$scale)]
  bn = data.frame(bn)
  mn = data.frame(rowMeans(bn, na.rm = T))
  colnames(mn) = clus[i]
  cm = cbind(cm, mn)
}

saveRDS(cm, 'cluster_mean_sex_effects_blocks.rds')
write.csv(cm, 'cluster_mean_sex_effects_blocks.csv') 

#### end ####

#### output significant genes per cluster ####

o = readRDS('blocks.rds')
clus = unique(o$sbm_cluster)
si = data.frame()
br = data.frame()
sr = sexres
sr$region = as.factor(sr$region)
levels(sr$region) = regions_short
sr$key = paste(sr$cell, sr$region, sep = ".")
for(i in 1:length(clus)){
  print(clus[i])
  bn = subset(sr, key %in% subset(o, sbm_cluster == clus[i])$scale)
  bn = data.frame(bn)
  sig = subset(bn, lfsr < 0.05)
  siglist <- sig %>%
    group_by(gene) %>%
    summarise(
      bias = if_else(
        all(bias == "male"), "male",
        if_else(all(bias == "female"), "female", "both")))
  siglist$cluster = clus[i]
  background = data.frame(gene = (unique(bn[complete.cases(bn),]$gene)))
  background$cluster = clus[i]
  si = rbind(si, siglist)
  br = rbind(br, background)
}

table(si$bias, si$cluster)
table(br$cluster)

write.csv(si, file = 'significant_genes_per_cluster.csv')
write.csv(br, file = 'background_genes_per_cluster.csv')

#### end ####

#### compare with sex-biased GMV ####

gmv = readRDS('sexGMV.rds')

# gmv sex bias similarity vs transcriptome sex bias similarity

sexc = readRDS('sexc.rds')

gmv_metrics <- c("bias", "absbias", "ranksigned", "rankabs")

all_names <- colnames(sexc)
celltypes <- unique(sub("\\..*", "", all_names))
regions <- gmv$region
region_pairs <- combn(regions, 2, simplify = FALSE)

compute_similarity_correlation <- function(celltype, metric) {

    full_names <- paste(celltype, regions, sep = ".")
  
  if (!all(full_names %in% colnames(sexc))) return(NULL)
  
  sim_matrix <- sexc[full_names, full_names]
  rownames(sim_matrix) <- colnames(sim_matrix) <- regions  
  
  sim_vals <- map_dbl(region_pairs, ~ sim_matrix[.x[1], .x[2]])
  gmv_vals <- map_dbl(region_pairs, ~ abs(
    gmv %>% filter(region == .x[1]) %>% pull(!!sym(metric)) -
      gmv %>% filter(region == .x[2]) %>% pull(!!sym(metric))
  ))
  
  # Correlation
  test <- cor.test(sim_vals, gmv_vals, method = "spearman")
  
  # Return
  tibble(
    celltype = celltype,
    metric = metric,
    correlation = test$estimate,
    p_value = test$p.value
  )
}

results <- expand_grid(celltype = celltypes, metric = gmv_metrics) %>%
  pmap_dfr(~ compute_similarity_correlation(..1, ..2))

results <- results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adj < 0.05 ~ "***",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

print(results)
write.csv(results, "gmv_pairwise_similarity_correlations.csv", row.names = FALSE)

results$celltype = factor(results$celltype, levels = c(exn.levels, inn.levels, other.levels))
ggplot(results, aes(x = celltype, y = correlation, color = metric)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.7)) +
  geom_text(aes(label = sig), vjust = -1, position = position_dodge(width = 0.7), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~metric, scales = "free_y") +
  scale_color_manual(values = c("bias" = "blue", "absbias" = "orange", "ranksigned" = "purple", "rankabs" = "green")) +
  labs(
    title = "Correlation of DEG proportions with Bias Metrics by Cell Type",
    x = "Cell Type",
    y = "Pearson Correlation (r)",
    color = "Metric") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

celltype <- "ADARB2"
metric <- "bias"

regions <- gmv$region
region_pairs <- combn(regions, 2, simplify = FALSE)
full_names <- paste(celltype, regions, sep = ".")

sim_matrix <- sexc[full_names, full_names]
rownames(sim_matrix) <- colnames(sim_matrix) <- regions

plot_data <- map_dfr(region_pairs, function(pair) {
  region1 <- pair[1]
  region2 <- pair[2]
  
  tibble(
    region1 = region1,
    region2 = region2,
    similarity = sim_matrix[region1, region2],
    gmv_diff = abs(
      gmv %>% filter(region == region1) %>% pull(!!sym(metric)) -
        gmv %>% filter(region == region2) %>% pull(!!sym(metric))
    )
  )
})

plot_data$pair = paste(plot_data$region1, plot_data$region2)
ggplot(plot_data, aes(x = gmv_diff, y = similarity)) +
  geom_point() + 
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_point(aes(x = gmv_diff, y = similarity, color = pair), size = 6) +
  geom_smooth(method='lm') +
  #scale_color_manual(values = region.colors) +
  theme_article()

# gmv vs mean deg props across celltypes

degs2 = degs %>% group_by(region) %>% summarise(prop = mean(prop))

res_joined <- degs2 %>%
  left_join(gmv, by = "region")

cor_results <- res_joined %>%
  summarise(
    cor_bias = cor(prop, bias, method = "spearman"),
    p_bias   = cor.test(prop, bias, method = "spearman")$p.value,
    
    cor_absbias = cor(prop, absbias, method = "spearman"),
    p_absbias   = cor.test(prop, absbias, method = "spearman")$p.value,
    
    cor_ranksigned = cor(prop, ranksigned, method = "spearman"),
    p_ranksigned   = cor.test(prop, ranksigned, method = "spearman")$p.value,
    
    cor_rankabs = cor(prop, rankabs, method = "spearman"),
    p_rankabs   = cor.test(prop, rankabs, method = "spearman")$p.value,
  )
View(cor_results)

ggplot(res_joined, aes(x = bias, y = prop)) +
  geom_point() + 
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_point(aes(x = bias, y = prop, color = region), size = 6) +
  geom_smooth(method='lm') +
  scale_color_manual(values = region.colors) +
  theme_article()

# gmv vs deg props within cell types

res_joined <- degs %>%
  left_join(gmv, by = "region")

cor_results <- res_joined %>%
  group_by(celltype) %>%
  summarise(
    cor_bias = cor(prop, bias, method = "spearman"),
    p_bias   = cor.test(prop, bias, method = "spearman")$p.value,
    
    cor_absbias = cor(prop, absbias, method = "spearman"),
    p_absbias   = cor.test(prop, absbias, method = "spearman")$p.value,
    
    cor_ranksigned = cor(prop, ranksigned, method = "spearman"),
    p_ranksigned   = cor.test(prop, ranksigned, method = "spearman")$p.value,
    
    cor_rankabs = cor(prop, rankabs, method = "spearman"),
    p_rankabs   = cor.test(prop, rankabs, method = "spearman")$p.value,
  )
View(cor_results)
write.csv(cor_results, file = 'gmv-prop-correlations.csv')

cor_results$p_bias_adj = p.adjust(cor_results$p_bias)
cor_results$p_absbias_adj = p.adjust(cor_results$p_absbias)
cor_results$p_ranksigned_adj = p.adjust(cor_results$p_ranksigned)
cor_results$p_rankabs_adj = p.adjust(cor_results$p_rankabs)

write.csv(cor_results, file = 'gmv-prop-correlations.csv')

library(tidyr)

cor_long <- cor_results %>%
  pivot_longer(
    cols = starts_with("cor_"),
    names_to = "metric",
    values_to = "correlation"
  )

pval_long <- cor_results %>%
  pivot_longer(
    cols = starts_with("p_"),    # this includes raw and adjusted p-values
    names_to = "metric_p",
    values_to = "p_value"
  ) %>%
  # Separate raw and adjusted p-values
  mutate(
    padj_flag = grepl("_adj$", metric_p),
    metric = sub("^p_", "cor_", metric_p),
    metric = sub("_adj$", "", metric)  # remove _adj suffix for matching
  ) %>%
  select(-metric_p)

pval_raw <- pval_long %>% filter(!padj_flag) %>% select(celltype, metric, p_value)
pval_adj <- pval_long %>% filter(padj_flag) %>% rename(p_value_adj = p_value) %>% select(celltype, metric, p_value_adj)

cor_pval_long <- cor_long %>%
  left_join(pval_raw, by = c("celltype", "metric")) %>%
  left_join(pval_adj, by = c("celltype", "metric")) %>%
  mutate(
    metric = recode(metric,
                    cor_bias = "bias",
                    cor_absbias = "absbias",
                    cor_ranksigned = "ranksigned",
                    cor_rankabs = "rankabs"),
    sig = case_when(
      p_value_adj < 0.05 ~ "***",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

ggplot(cor_pval_long, aes(x = celltype, y = correlation, color = metric)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.7)) +
  geom_text(aes(label = sig), position = position_dodge(width = 0.7), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~metric, scales = "free_y") +
  scale_color_manual(values = c("bias" = "blue", "absbias" = "orange", "ranksigned" = "purple", "rankabs" = "green")) +
  labs(
    title = "Correlation of DEG proportions with Bias Metrics by Cell Type",
    x = "Cell Type",
    y = "Pearson Correlation (r)",
    color = "Metric") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# plot example

cellnow = 'SST'
measure = 'bias'

pl = subset(res_joined, celltype == cellnow)[,c('region','prop', measure)]
colnames(pl)[3] = 'measure'
ggplot(pl, aes(x = measure, y = prop)) +
  geom_point() + 
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_point(aes(x = measure, y = prop, color = region), size = 6) +
  geom_smooth(method='lm') +
  scale_color_manual(values = region.colors) +
  theme_article()

#### end ####

#### plot mashr sex effects (chr enrich) ####

hsap = useEnsembl(biomart = 'ensembl', dataset='hsapiens_gene_ensembl', version = 98) 
gene2chromall = getBM(attributes=c('external_gene_name','chromosome_name','gene_biotype'), mart = hsap)
gene2chromall = subset(gene2chromall, chromosome_name %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
gene2chromall = subset(gene2chromall, 
                       gene_biotype %in% c('protein_coding','lncRNA',
                                           'IG_C_gene','IG_D_gene','IG_J_gene',
                                           'IG_LV_gene','IG_V_gene','IG_V_pseudogene',
                                           'IG_J_pseudogene','IG_C_pseudogene','TR_C_gene',
                                           'TR_D_gene','TR_J_gene','TR_V_gene',
                                           'TR_V_pseudogene','TR_J_pseudogene'))

colnames(gene2chromall)[1] = 'gene'
length(unique(gene2chromall$gene))

pl = merge(sexres, gene2chromall, by = 'gene', all = T)
pl = pl[complete.cases(pl),]
length(unique(pl$gene))

result_df = pl %>%
  group_by(gene, chromosome_name) %>%
  summarise(
    label = ifelse(all(bias == 'ns'), 'ns', 'bias'))
table(result_df$chromosome_name)

proportion_df <- result_df %>%
  group_by(chromosome_name, label) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(chromosome_name) %>%
  mutate(proportion = count / sum(count))

proportion_df$chromosome_name = factor(proportion_df$chromosome_name,
                                      levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
ggplot(proportion_df, aes(x = chromosome_name, y = proportion, fill = label)) +
  geom_bar(stat = 'identity', position = 'stack')

result_df = pl %>%
  group_by(gene, chromosome_name) %>%
  summarise(
    label = case_when(
      all(bias == 'ns') ~ 'ns',
      all(bias %in% c('male', 'ns')) ~ 'male',
      all(bias %in% c('female', 'ns')) ~ 'female',
      any(bias == 'female') & any(bias == 'male') ~ 'both'
    )
  )

proportion_df <- result_df %>%
  group_by(chromosome_name, label) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(chromosome_name) %>%
  mutate(proportion = count / sum(count))

proportion_df$chromosome_name = factor(proportion_df$chromosome_name,
                                       levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
proportion_df$label = factor(proportion_df$label, levels = c('female','male','both','ns'))
ggplot(proportion_df, aes(x = chromosome_name, y = proportion, fill = label)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c(sex.colors,'darkgreen','#666666')) +
  theme_article()

#### end ####

#### plot mashr sex effects (x chromosome) ####

xci = read.csv('xci-status.csv')
xci$Gene.ID = str_sub(xci$Gene.ID,1,15)
# missing TMSB4X and SOX3
# gene names have changed

colnames(xci)[1] = 'gene'
colnames(xci)[2] = 'ensembl_gene_id'
pl = sexres[complete.cases(sexres),]
check = merge(gene2chrom, xci, by = 'ensembl_gene_id', all = T)
check = subset(check, chromosome_name == 'X')
View(subset(check, gene.x != gene.y & gene.x %in% pl$gene))
toadd = subset(check, gene.x != gene.y & gene.x %in% pl$gene)
add = data.frame(gene = toadd$gene.x, ensembl_gene_id = toadd$ensembl_gene_id, 
                 Transcript.type='protein_coding', Combined.XCI.status = toadd$Combined.XCI.status)
xci = rbind(xci, add)
add = data.frame(gene = 'TMSB4X',ensembl_gene_id='ENSG00000205542',Transcript.type='protein_coding',Combined.XCI.status='escape')
xci = rbind(xci, add)
add = data.frame(gene = 'SOX3',ensembl_gene_id='ENSG00000134595',Transcript.type='protein_coding',Combined.XCI.status='escape')
xci = rbind(xci, add)

pl = subset(pl, gene %in% xci$gene)
pl$gene = droplevels(pl$gene)
pl = pl %>% complete(gene, cell, region)
pl = merge(pl, xci, by = 'gene', all.x = T)
#pl = subset(pl, lfsr < 0.05)
pl$beta = ifelse(pl$beta < -0.5, -0.5, pl$beta)
pl$beta = ifelse(pl$beta > 0.5, 0.5, pl$beta)
pl$Combined.XCI.status = ifelse(pl$gene %in% par$external_gene_name, 'PAR', pl$Combined.XCI.status)
pl$region = factor(pl$region)
levels(pl$region) = regions_short
pl$cell = factor(pl$cell, levels = c(exn.levels, inn.levels, other.levels))
pl$beta = ifelse(pl$lfsr > 0.05, 0, pl$beta)

ggplot(pl, aes(x = gene, y = cell, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(region~Combined.XCI.status, scales = 'free', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 2))

pl = subset(pl, Combined.XCI.status != 'inactive')

ggplot(pl, aes(x = gene, y = cell, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(region~Combined.XCI.status, scales = 'free', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 2))

gams = c('PRKX','PRKY','ZFX','ZFY','USP9X','USP9Y','TBL1X','TBL1Y','NLGN4Y','NLGN4X','UTY','KDM6A','TMSB4X','TMSB4Y','TGIF2LX','TGIF2LY',
         'RPS4X','RPS4Y1','DDX3Y','DDX3X','EIF1AY','EIF1AX','KDM5C','KDM5D','PCDH11Y','PCDH11X','TXLNG','TXLNGY','AMELX','AMELY','SOX3','SRY')

pl = subset(pl, Combined.XCI.status == 'escape' | Combined.XCI.status == 'PAR' | gene %in% gams)
pl$gam = ifelse(pl$gene %in% gams, 'gam', pl$Combined.XCI.status)
pl$gam = factor(pl$gam, levels= c('gam','escape','PAR'))
levels(pl$gam)

ggplot(pl, aes(x = cell, y = gene, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(gam~region, scales = 'free', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        strip.text = element_text(size = 18),
        axis.title = element_blank())

pl = sexres[complete.cases(sexres),]
pl = merge(pl, gene2chrom[,c('gene','chromosome_name')], by = 'gene', all.x = T)
pl = merge(pl, xci[,c('gene','Combined.XCI.status')], by = 'gene', all.x = T)
table(pl$chromosome_name, pl$Combined.XCI.status)
pl$key = ifelse(pl$chromosome_name == 'X', paste("X",pl$Combined.XCI.status,sep="-"), pl$chromosome_name)
table(pl$key)
pl = subset(pl, lfsr < 0.05)
pl$beta = ifelse(pl$beta < -5, -5, pl$beta)
pl$beta = ifelse(pl$beta > 5, 5, pl$beta)
pl$key = factor(pl$key, levels = c('Y','PAR','X-escape','X-variable','X-inactive',22:1))
pl = pl[complete.cases(pl$key),]

ggplot() +
  geom_violin(data = subset(pl, beta > 0), aes(x = beta, y = key), fill = sex.colors[2], width = 1) +
  geom_violin(data = subset(pl, beta < 0), aes(x = beta, y = key), fill = sex.colors[1], width = 1) +
  geom_hline(yintercept = 0, linetype ='dashed') +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

pl$key = factor(pl$key, levels = c(1:22,'X-inactive','X-escape','PAR','Y'))

ggplot() +
  geom_violin(data = subset(pl, beta > 0), aes(y = beta, x = key), fill = sex.colors[2], width = 1) +
  geom_violin(data = subset(pl, beta < 0), aes(y = beta, x = key), fill = sex.colors[1], width = 1) +
  geom_vline(xintercept = 0, linetype ='dashed') +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

# balaton

xci = read.csv('xci-balaton.csv')
colnames(xci)[1] = 'gene'

pl = sexres[complete.cases(sexres),]
pl = subset(pl, gene %in% xci$gene)
pl$gene = droplevels(pl$gene)
pl = pl %>% complete(gene, cell, region)
pl = merge(pl, xci, by = 'gene', all.x = T)
#pl = subset(pl, lfsr < 0.05)
pl$beta = ifelse(pl$beta < -0.5, -0.5, pl$beta)
pl$beta = ifelse(pl$beta > 0.5, 0.5, pl$beta)
pl$region = factor(pl$region)
levels(pl$region) = regions_short
pl$cell = factor(pl$cell, levels = c(exn.levels, inn.levels, other.levels))
pl$beta = ifelse(pl$lfsr > 0.05, 0, pl$beta)

ggplot(pl, aes(x = cell, y = gene, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(Domain.category~region, scales = 'free', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=6))

gams = c('PRKX','PRKY','ZFX','ZFY','USP9X','USP9Y','TBL1X','TBL1Y','NLGN4Y','NLGN4X','UTY','KDM6A','TMSB4X','TMSB4Y','TGIF2LX','TGIF2LY',
         'RPS4X','RPS4Y1','DDX3Y','DDX3X','EIF1AY','EIF1AX','KDM5C','KDM5D','PCDH11Y','PCDH11X','TXLNG','TXLNGY','AMELX','AMELY','SOX3','SRY')

pl$gam = ifelse(pl$gene %in% gams, 'gam', pl$Domain.category)

ggplot(pl, aes(x = cell, y = gene, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(gam~region, scales = 'free', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.y.right = element_text(angle = 0, size = 10,  hjust = 1),
        axis.title = element_blank())

pl = sexres[complete.cases(sexres),]
pl = merge(pl, gene2chrom[,c('gene','chromosome_name')], by = 'gene', all.x = T)
pl = merge(pl, xci[,c('gene','Domain.category')], by = 'gene', all.x = T)
table(pl$chromosome_name, pl$Domain.category)
pl$key = ifelse(pl$chromosome_name == 'X', paste("X",pl$Combined.XCI.status,sep="-"), pl$chromosome_name)
table(pl$key)
pl = subset(pl, lfsr < 0.05)
pl$beta = ifelse(pl$beta < -5, -5, pl$beta)
pl$beta = ifelse(pl$beta > 5, 5, pl$beta)
pl$key = factor(pl$key, levels = c('Y','PAR','X-escape','X-variable','X-inactive',22:1))
pl = pl[complete.cases(pl$key),]

ggplot() +
  geom_violin(data = subset(pl, beta > 0), aes(x = beta, y = key), fill = sex.colors[2], width = 1) +
  geom_violin(data = subset(pl, beta < 0), aes(x = beta, y = key), fill = sex.colors[1], width = 1) +
  geom_hline(yintercept = 0, linetype ='dashed') +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

pl$key = factor(pl$key, levels = c(1:22,'X-inactive','X-escape','PAR','Y'))

ggplot() +
  geom_violin(data = subset(pl, beta > 0), aes(y = beta, x = key), fill = sex.colors[2], width = 1) +
  geom_violin(data = subset(pl, beta < 0), aes(y = beta, x = key), fill = sex.colors[1], width = 1) +
  geom_vline(xintercept = 0, linetype ='dashed') +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

# all X chr genes

xgenes = subset(gene2chrom, chromosome_name == 'X')

pl = sexres[complete.cases(sexres),]
pl = subset(pl, gene %in% xgenes$gene)
pl$gene = droplevels(pl$gene)
pl = pl %>% complete(gene, cell, region)
pl = merge(pl, xci, by = 'gene', all.x = T)
pl$beta = ifelse(pl$beta < -0.5, -0.5, pl$beta)
pl$beta = ifelse(pl$beta > 0.5, 0.5, pl$beta)
pl$Combined.XCI.status = ifelse(pl$gene %in% par$external_gene_name, 'PAR', pl$Combined.XCI.status)
pl$Combined.XCI.status[is.na(pl$Combined.XCI.status)] = 'nonXCI'
pl$region = factor(pl$region)
levels(pl$region) = regions_short
pl$cell = factor(pl$cell, levels = c(exn.levels, inn.levels, other.levels))
pl$beta = ifelse(pl$lfsr > 0.05, 0, pl$beta)

table(pl$Combined.XCI.status)

gams = c('PRKX','PRKY','ZFX','ZFY','USP9X','USP9Y','TBL1X','TBL1Y','NLGN4Y','NLGN4X','UTY','KDM6A','TMSB4X','TMSB4Y','TGIF2LX','TGIF2LY',
         'RPS4X','RPS4Y1','DDX3Y','DDX3X','EIF1AY','EIF1AX','KDM5C','KDM5D','PCDH11Y','PCDH11X','TXLNG','TXLNGY','AMELX','AMELY','SOX3','SRY')
pl = subset(pl, Combined.XCI.status == 'escape' | Combined.XCI.status == 'PAR' | gene %in% gams)
pl$gam = ifelse(pl$gene %in% gams, 'gam', pl$Combined.XCI.status)

#### end ####

#### rank DEGs vs sexTWI ####

# DEGs

sex = readRDS('/mashr/mashr_sex_class_approach3.rds')
sexbeta = get_pm(sex)
sexlfsr = get_lfsr(sex)
sexbeta2 = reshape2::melt(sexbeta)
sexlfsr2 = reshape2::melt(sexlfsr)
sexres = cbind(sexbeta2, sexlfsr2)
d = str_split_fixed(sexres$Var2, pattern = '\\.', 2)
sexres = cbind(sexres, d)
sexres = sexres[,c(1,3,6:8)]
colnames(sexres) = c('gene','beta','lfsr','cell','region')
sexres$key = paste(sexres$region, sexres$cell, sexres$gene)
#sexres$beta = ifelse(sexres$key %in% diffsign$key, 0, sexres$beta)
#sexres$lfsr = ifelse(sexres$key %in% diffsign$key, 1, sexres$lfsr)
sexres$bias = ifelse(sexres$beta > 0 & sexres$lfsr < 0.05, 'male', 'ns')
sexres$bias = ifelse(sexres$beta < 0 & sexres$lfsr < 0.05, 'female', sexres$bias)
sexres$key = paste(sexres$region, sexres$cell)

meta = readRDS('meta_merged.rds')
meta$sample_id2 = as.character(meta$sample_id2)
check = reshape2::melt(meta[,c(3,12,44:154)]) 
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var2, check$Var2)
remove = subset(check, Freq < 10)
sexres = subset(sexres, key %!in% remove$key)

pl = sexres[complete.cases(sexres),]
pl = merge(pl, unique(gene2chrom[,c('gene','chromosome_name')]), by = 'gene', all.x = T)
unique(pl[which(is.na(pl$chromosome_name)),'gene'])
pl = pl[complete.cases(pl),]
pl$chrom2 = ifelse(pl$chromosome_name == 'X' | pl$chromosome_name == 'Y' | pl$chromosome_name == 'PAR', pl$chromosome_name, 'Autosomal')
pl = pl %>% group_by(cell, region, bias, chrom2) %>% summarise(n = n())
pl$key = paste(pl$cell, pl$region)
pl2 = pl %>% group_by(cell, region) %>% summarise(sum = sum(n))
pl2$key = paste(pl2$cell, pl2$region)
pl = merge(pl, pl2, by = 'key')
pl$prop = pl$n / pl$sum
pl$cell.x = factor(pl$cell.x,levels = c(exn.levels, inn.levels, other.levels))
pl$region.x = factor(pl$region.x)
levels(pl$region.x) = regions_short

pl = subset(pl, bias != 'ns')

degs = pl %>% group_by(cell.x, region.x) %>% summarise(prop = sum(prop), n = sum(n))
colnames(degs) = c('celltype','region','prop','n')
degs$rankdeg = rank(degs$n, ties.method = "random")  
degs$rankprop = rank(degs$prop, ties.method = "random")  

mod = aov(prop ~ celltype + region, degs)
summary(mod)
TukeyHSD(mod)

# twi

res = readRDS('/TRADE/all_SEX_TRADE_class_outliers.rds')
colnames(res) = c('twi','region','celltype')

meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,44:154)]) 
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
res$key = paste(res$region, res$celltype)
res = subset(res, key %!in% remove$key)
res$region = factor(res$region)
levels(res$region) = regions_short
res$ranktwi = rank(res$twi, ties.method = "random")  

pl = merge(degs, res, by = c('celltype','region'))

ggplot(pl, aes(x = rankdeg, y = rankprop)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  theme_article()

ggplot(pl, aes(x = rankdeg, y = ranktwi)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  theme_article()

ggplot(pl, aes(x = rankprop, y = ranktwi)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  theme_article()

cor.test(pl$rankdeg, pl$rankprop)
cor.test(pl$rankdeg, pl$ranktwi)
cor.test(pl$rankprop, pl$ranktwi)

#### end ####


