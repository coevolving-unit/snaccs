library(speckle)
library(limma)
library(ggplot2)
library(statmod)

#### load data ####

meta = readRDS('meta_merged.rds')
meta = subset(meta, bad == '#N/A')
meta = subset(meta, repeat. != 'repeat')
dim(meta)

#### end ####

#### Figure 1C ####

## plot counts subclasses x regions

pl = data.frame()
for(i in 1:length(regions)){
  m = subset(meta, region == regions[i])
  pln = data.frame(total = colSums(m[,c(44:67)]), 
                   cell = colnames(m)[c(44:67)])
  pln$cell = factor(pln$cell,levels = c(exn.levels, inn.levels, other.levels))
  pln$region = regions[i]
  pl = rbind(pl, pln)
}
levels(pl$cell)
pl$region = factor(pl$region)
levels(pl$region) = regions_short

ggplot(pl, aes(x = cell, y = total, fill = cell)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(exn.colors, inn.colors, other.colors)) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

#### end ####

#### plot sex diffs raw proportions ####

pl = data.frame()
for(i in 1:length(regions)){
  mmeta = subset(meta, sex == 'Male' & region == regions[i])
  fmeta = subset(meta, sex == 'Female' & region == regions[i])
  m = data.frame(total = colSums(mmeta[,c(44:67)]), cell = colnames(mmeta)[c(41:64)])
  f = data.frame(total = colSums(fmeta[,c(44:67)]), cell = colnames(fmeta)[c(41:64)])
  pln = rbind(cbind(m, sex = 'm'), cbind(f,sex='f'))
  pln$cell = factor(pln$cell,
                    levels = c(exn.levels, inn.levels, other.levels))
  pln$region = regions[i]
  pl = rbind(pl, pln)
}
levels(pl$cell)
pl$region = factor(pl$region)
levels(pl$region) = regions_short

ggplot(pl, aes(x = cell, y = total, fill = sex)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = sex.colors) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

meta2 = reshape2::melt(meta[,c(1,3,19,7:11,20,42,43,68,44:67)], 
                       id.vars = c('new.id','region','RIN','PMI','ph','individual','sex','age','batch','nreads','maprate','total_cells'))

metap = meta2 %>% group_by(new.id) %>% mutate(freq = value / sum(value))
metap$variable = factor(metap$variable,
                        levels = c(exn.levels, inn.levels, other.levels))
pl = metap %>% group_by(region, sex, variable) %>% summarise(mean = mean(freq))
pl$region = factor(pl$region)
levels(pl$region) = regions_short

ggplot(pl, aes(x = variable, y = mean, fill = sex)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = sex.colors) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

ggplot(pl, aes(x = variable, y = mean, fill = sex)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = sex.colors) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

#### end ####

#### calc sex diffs within each region ####

meta2 = reshape2::melt(meta[,c(1,3,19,7:11,20,42,43,68,44:67)], 
                       id.vars = c('new.id','region','RIN','PMI','ph','individual','sex','age','batch','nreads','maprate','total_cells'))

cells = data.frame(celltype = rep(meta2$variable, meta2$value))
cells$new.id = rep(meta2$new.id, meta2$value)
cells$region = rep(meta2$region, meta2$value)
cells$individual = rep(meta2$individual, meta2$value)
cells$batch = rep(meta2$batch, meta2$value)
cells$sex = rep(meta2$sex, meta2$value)
cells$age = rep(meta2$age, meta2$value)
cells$PMI = rep(meta2$PMI, meta2$value)
cells$RIN = rep(meta2$RIN, meta2$value)
cells$nreads = rep(meta2$nreads, meta2$value)
cells$total_cells = rep(meta2$total_cells, meta2$value)
cells$key = paste(cells$region, cells$celltype)
dim(cells)

table(cells$celltype)

out = data.frame()
out2 = data.frame()
alt_prop_all = data.frame()

for(i in 1:length(regions)){
  
  # unadjusted 
  # proportions ~ sex
  
  datnow = subset(cells, region == regions[i])
  pro = propeller(x = datnow, clusters = datnow$celltype, sample = datnow$new.id, group = datnow$sex)
  rownames(pro) = NULL
  pro$region = regions_short[i]
  out = rbind(out, pro)
  
  # adjusted for covariates
  # proportion ~ sex + age + PMI + RIN + batch + nreads + total cells
  
  alt_prop = getTransformedProps(clusters = datnow$celltype, sample = datnow$new.id)
  alt_prop_all = rbind(alt_prop_all, alt_prop$TransformedProps)
  sample_data_now = subset(cells, new.id %in% datnow$new.id)
  sample_data_now = sample_data_now[,c(2:11)]
  sample_data_now = unique(sample_data_now)
  idx = match(colnames(alt_prop$TransformedProps), sample_data_now$new.id)
  sample_data_now = sample_data_now[idx,]
  design = model.matrix(~ 0 + sample_data_now$sex + 
                          sample_data_now$age + 
                          sample_data_now$PMI +
                          sample_data_now$RIN + 
                          sample_data_now$batch + 
                          sample_data_now$nreads + 
                          sample_data_now$total_cells)
  cont.sex <- c(1,-1,rep(0,dim(design)[2]-2)) # female vs male
  pro2 = propeller.ttest(alt_prop, design = design, contrasts=cont.sex,
                         trend = FALSE, robust = TRUE, sort=TRUE)
  pro2$celltype = rownames(pro2)
  rownames(pro2) = NULL
  pro2$region = regions_short[i]
  out2 = rbind(out2, pro2)
}

View(out)
View(out2)

saveRDS(out, file = '/proportions/proportions-within-region-unadjusted.rds')
saveRDS(out2, file = '/proportions/proportions-within-region-adjusted.rds')

write.csv(out2, file = '/proportions/proportions-within-region-adjusted.csv')

# glia vs neurons

cells$group = ifelse(cells$celltype %in% other.levels, 'glia', 'neuron')
table(cells$group)

out = data.frame()
out2 = data.frame()
alt_prop_all = data.frame()

for(i in 1:length(regions)){
  
  # unadjusted 
  # proportions ~ sex
  
  datnow = subset(cells, region == regions[i])
  pro = propeller(x = datnow, clusters = datnow$group, sample = datnow$new.id, group = datnow$sex)
  rownames(pro) = NULL
  pro$region = regions_short[i]
  out = rbind(out, pro)
  
  # adjusted for covariates
  # proportion ~ sex + age + PMI + RIN + batch + nreads + total cells
  
  alt_prop = getTransformedProps(clusters = datnow$group, sample = datnow$new.id)
  alt_prop_all = rbind(alt_prop_all, alt_prop$TransformedProps)
  sample_data_now = subset(cells, new.id %in% datnow$new.id)
  sample_data_now = unique(sample_data_now)
  idx = match(colnames(alt_prop$TransformedProps), sample_data_now$new.id)
  sample_data_now = sample_data_now[idx,]
  design = model.matrix(~ 0 + sample_data_now$sex + 
                          sample_data_now$age + 
                          sample_data_now$PMI +
                          sample_data_now$RIN + 
                          sample_data_now$batch + 
                          sample_data_now$nreads + 
                          sample_data_now$total_cells)
  cont.sex <- c(1,-1,rep(0,dim(design)[2]-2)) # female vs male
  pro2 = propeller.ttest(alt_prop, design = design, contrasts=cont.sex,
                         trend = FALSE, robust = TRUE, sort=TRUE)
  pro2$celltype = rownames(pro2)
  rownames(pro2) = NULL
  pro2$region = regions_short[i]
  out2 = rbind(out2, pro2)
}

View(out)
View(out2)

#### end ####

#### calc sex diffs across regions ####

# individual included as a random effect

pro = propeller(x = cells, clusters = cells$celltype, sample = cells$new.id, group = cells$sex)
#pro = propeller(x = cells, clusters = cells$group, sample = cells$new.id, group = cells$sex)

rownames(pro) = NULL
View(pro)

saveRDS(pro, file = '/proportions/proportions-across-regions-unadjusted.rds')

alt_prop = getTransformedProps(clusters = cells$celltype, sample = cells$new.id)
#alt_prop = getTransformedProps(clusters = cells$group, sample = cells$new.id)

sample_data_now = cells
sample_data_now = sample_data_now[,c(2:11)]
sample_data_now = unique(sample_data_now)
idx = match(colnames(alt_prop$TransformedProps), sample_data_now$new.id)
sample_data_now = sample_data_now[idx,]

design = model.matrix(~ 0 + sample_data_now$sex + 
                        sample_data_now$age + 
                        sample_data_now$PMI +
                        sample_data_now$RIN + 
                        sample_data_now$batch + 
                        sample_data_now$nreads + 
                        sample_data_now$total_cells+
                        sample_data_now$region)
colnames(design) = gsub("sample_data_now","",colnames(design))
colnames(design) = gsub("\\$","",colnames(design))

dupcor = duplicateCorrelation(alt_prop$TransformedProps, design=design, block=sample_data_now$individual)
fit1 = lmFit(alt_prop$TransformedProps, design=design, block=sample_data_now$individual, correlation=dupcor$consensus)
contr = makeContrasts(sexMale - sexFemale, levels = colnames(coef(fit1)))
tmp = contrasts.fit(fit1, contr)
tmp = eBayes(tmp)
pro2 = data.frame(tmp$coefficients)
pro2$p = c(tmp$p.value)
pro2$celltype = rownames(pro2)
rownames(pro2) = NULL
colnames(pro2) = c('beta','p','celltype')
pro2$FDR = p.adjust(pro2$p, method = 'fdr')
View(pro2)

saveRDS(pro2, file = '/proportions/proportions-across-regions-adjusted.rds')

write.csv(pro2, file = '/proportions/proportions-across-regions-adjusted.csv')

# proportion ~ cell type * sex' with random term (~1|person/region)

design = model.matrix(~ 0 + cells$sex * cells$celltype)
dupcor = duplicateCorrelation(alt_prop$TransformedProps, design=design, block=cells$individual)
fit1 = lmFit(alt_prop$TransformedProps, design=design, block=sample_data_now$individual, correlation=dupcor$consensus)
contr = makeContrasts(sexMale - sexFemale, levels = colnames(coef(fit1)))
tmp = contrasts.fit(fit1, contr)
tmp = eBayes(tmp)

# Figure 2A - plot adjusted proportions

p = apply(as.matrix.noquote(t(alt_prop$TransformedProps)),2,as.numeric)
rownames(p) = colnames(alt_prop$TransformedProps)
dim(p)
dim(sample_data_now)

res = removeBatchEffect(t(p), 
                  batch = sample_data_now$batch, 
                  batch2 = sample_data_now$region,
                  covariates = sample_data_now[,c('age','PMI','RIN','nreads','total_cells')])

res = melt(res)
colnames(res) = c('celltype','new.id','prop')
res = merge(res, unique(meta[,c('new.id','sex','region')]), all.x = T)
res$celltype = factor(res$celltype, levels = c(exn.levels, inn.levels, other.levels))
res$region = factor(res$region)
levels(res$region) = regions_short

res$key = paste(res$region, res$celltype)
meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,41:151)])
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$Var1 = as.factor(check$Var1)
levels(check$Var1) = regions_short
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
res = subset(res, key %!in% remove$key)

ggplot(res, aes(x = celltype, y = prop, fill = sex)) +
  geom_violin(scale = "width") +
  geom_boxplot(outlier.colour = NA, ) +
  scale_fill_manual(values = sex.colors) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_blank())

pl = res %>% group_by(region, sex, celltype) %>% summarise(mean = mean(prop))
pl = pl[complete.cases(pl),]

ggplot(pl, aes(x = celltype, y = mean, fill = sex)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = sex.colors) +
  facet_wrap(~region, ncol = 1) +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())

# Figure 2A - plot fold changes

#pl = readRDS('/proportions/proportions-across-regions-adjusted.rds')
pl = readRDS('/proportions/proportions-within-region-adjusted.rds')

pl$key = paste(pl$region, pl$celltype)
meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,41:151)])
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$Var1 = as.factor(check$Var1)
levels(check$Var1) = regions_short
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
pl = subset(pl, key %!in% remove$key)

pl$log2 = log2(pl$PropRatio)
pl$col = ifelse(pl$log2 < 0, 'female', 'male')
pl$celltype = factor(pl$celltype, levels = c(exn.levels, inn.levels, other.levels))

ggplot(pl, aes(x = celltype, y = log2, fill = col)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = sex.colors) +
  facet_wrap(~region, ncol = 1, drop = TRUE) +
  theme_article() +
  ylim(c(-1.5, 1.5)) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_blank())


#### end ####

