library(variancePartition)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(corrplot)
library(reshape2)
library(egg)

#### overall variance partitioning  ####

v = readRDS('/variance_partitioning/varpart_all_samples_class_corrected_outliers.rds')

write.csv(v, file = 'overall_variance_partitioning.csv')

vc = cor(v, method = 'spearman')
vc = vc[1:9,1:9]
rval = psych::corr.test(v[,c(1:9)], method = 'spearman', adjust="none")$r
write.csv(rval, 'overall-var-part-cor-rho.csv')
pval = psych::corr.test(v[,c(1:9)], method = 'spearman', adjust="none")$p
write.csv(pval, 'overall-var-part-cor-p.csv')
pval = psych::corr.test(v[,c(1:9)], method = 'spearman', adjust="fdr")$p
write.csv(pval, 'overall-var-part-cor-fdr.csv')
colnames(vc) = rownames(vc) = colnames(pval) = rownames(pval) = 
          c('cell type','individual','region','sex','age','RIN','PMI','# reads','# cells')

corrplot(vc, type="upper", order = 'AOE',tl.cex = 2, tl.pos="l", tl.col = "black")
corrplot(vc, type="lower", p.mat=pval, insig='p-value', 
         sig.level=0, add=T, order = 'AOE', 
         tl.pos="t", tl.cex = 2, tl.col = "black", cl.pos="n",
         col = COL2('PiYG'))

corrplot(vc, order = 'AOE', 
         tl.pos="t", tl.cex = 2, tl.col = "black", 
         col = COL2('PRGn'))

v$gene = rownames(v)
v = melt(v)
v$variable = factor(v$variable, levels = c('celltype','individual','total_cells',
                                           'region','sex','PMI','RIN','age','nreads','Residuals'))
levels(v$variable) = c('cell type','individual','# cells','region',
                       'sex','PMI','RIN','age','# reads','residuals')

ggplot(v, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_article() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gene2chrom = readRDS('gene2chrom_mash.rds')
table(gene2chrom$chromosome_name)
colnames(gene2chrom)[2] = 'gene'
gene2chrom = unique(gene2chrom[,c('gene','chromosome_name')])
gene2chrom = subset(gene2chrom, chromosome_name %in% c(1:22, 'X', 'Y'))

par = readRDS('par.rds')
gene2chrom[which(gene2chrom$gene %in% par$external_gene_name),'chromosome_name'] = 'PAR'

vp.chrom = v
vp.chrom = merge(vp.chrom, gene2chrom, by = 'gene')
vp.chrom$chrom2 = ifelse(vp.chrom$chrom == 'X' | vp.chrom$chrom == 'Y' | vp.chrom$chrom == 'PAR', vp.chrom$chrom, 'Autosomal')
data.frame(vp.chrom %>% group_by(variable, chrom2) %>% summarise(mean = mean(value)*100))
mean = data.frame(v %>% group_by(variable) %>% summarise(mean = mean(value)))
mean$mean = round(mean$mean, digits = 3)

g = ggplot(v, aes(x=variable, y=value)) + 
  geom_boxplot()
split <- split(vp.chrom, vp.chrom$variable)
ld <- layer_data(g)
outliers <- lapply(seq_along(split), function(i) {
  box <- ld[ld$group == i, ]
  data <- split[[i]]
  data <- data[data$value > box$ymax , ]
  data
})
outliers <- do.call(rbind, outliers)

getPalette = colorRampPalette(brewer.pal(11, "Dark2"))
colnow = getPalette(12)

ggplot(v, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(outlier.alpha = 0.1, alpha = 0.6) + 
  scale_fill_manual(values = colnow) +
  theme_article() +
  geom_text(data = mean, aes(label = mean, y = -0.02), size = 4) +
  ylab('proportion of variance explained') +
  geom_point(data = outliers, aes(shape = chrom2, size = chrom2)) +
  scale_shape_manual(values = c(1, 15, 4, 17)) +
  scale_size_manual(values = c(NA, 2, 2, 2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = 'none',
        axis.title.x = element_blank())

chr = vp.chrom %>% group_by(variable, chrom2) %>% summarise(mean = mean(value)) 
View(chr)

#### end  ####

#### per cell subclass ####

v = readRDS('/variance_partitioning/all_varpart_continuous_corrected_class.rds')

ev = subset(v, celltype %in% exn.levels)
ev[,c(1:7)] = round(ev[,c(1:7)], 3)
iv = subset(v, celltype %in% inn.levels)
iv[,c(1:7)] = round(iv[,c(1:7)], 3)
ov = subset(v, celltype %in% other.levels)
ov[,c(1:7)] = round(ov[,c(1:7)], 3)
write.csv(ev, file = 'exn_variance_partitioning_continuous.csv')
write.csv(iv, file = 'inn_variance_partitioning_continuous.csv')
write.csv(ov, file = 'other_variance_partitioning_continuous.csv')

v$celltype = factor(v$celltype, levels = c(exn.levels, inn.levels, other.levels))
vp.chrom = v
vp.chrom = merge(vp.chrom, gene2chrom, by = 'gene')
vp.chrom$chrom2 = ifelse(vp.chrom$chrom == 'X' | vp.chrom$chrom == 'Y' | vp.chrom$chrom == 'PAR', vp.chrom$chrom, 'Autosomal')

## filter

meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,44:154)]) 
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
vp.chrom$key = paste(vp.chrom$region, vp.chrom$celltype)
vp.chrom = subset(vp.chrom, key %!in% remove$key)
vp.chrom$region = as.factor(vp.chrom$region)
levels(vp.chrom$region) = regions_short

obj = data.frame(vp.chrom)
g <- ggplot(obj, aes(x=sex, y=celltype)) + geom_boxplot() + facet_wrap(~region, ncol = 1)
obj$key = paste(obj$region, obj$celltype, sep = ".")
split <- split(obj, list(obj$region, obj$celltype))
split = split[which(names(split) %in% obj$key)]
ld <- layer_data(g)
ld$code = c(1:dim(ld)[1])
outliers <- lapply(seq_along(split), function(i) {
  box <- ld[ld$code == i, ]
  data <- split[[i]]
  data <- data[data$sex > box$xmax , ]
  data
})
outliers <- do.call(rbind, outliers)

vp.chrom <- vp.chrom %>%
  mutate(celltype_region = interaction(region, celltype)) %>%
  group_by(celltype_region) %>%
  mutate(median_sex = median(sex, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(celltype_region = fct_reorder(celltype_region, median_sex))

library(tidytext)  # for reorder_within and scale_x_reordered

ggplot(vp.chrom, aes(
  x = reorder_within(celltype, -sex, region, FUN = median),
  y = sex,
  fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.1, alpha = 0.6) + 
  facet_wrap(~region, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = c(exn.colors, inn.colors, other.colors)) +
  theme_article() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_blank()) +
  geom_point(data = outliers, aes(
    x = reorder_within(celltype, sex, region, FUN = median),
    y = sex,
    shape = chrom2,
    size = chrom2)) +
  scale_shape_manual(values = c(1, 15, 4, 17)) +
  scale_size_manual(values = c(NA, 1, 1, 1)) +
  scale_x_reordered()

p = vp.chrom %>% group_by(celltype) %>% summarise(mean = mean(sex))
View(p)

p = vp.chrom %>% group_by(region) %>% summarise(mean = mean(sex))
View(p)

p = vp.chrom %>% group_by(celltype, region) %>% summarise(mean = mean(sex))
mean(p$mean)
min(p$mean)
max(p$mean)
View(p)

#### end ####

#### correlations ####

celltypes = unique(vp.chrom$celltype)
out = data.frame()
for(i in 1:length(regions_short)){
  for(j in 1:length(celltypes)){
    now = subset(vp.chrom, region == regions_short[i] & celltype == celltypes[j])
    c = cor(now[,c(2:8)])
    c = melt(c)
    c$region = regions_short[i] 
    c $celltype = celltypes[j]
    out= rbind(out, c)
}}
out = subset(out, value != 1)

write.csv(out, file = 'var-part-correlations-within-regions-celltypes.csv')

#### end ####
