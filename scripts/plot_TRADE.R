library(ggplot2)
library(egg)
library(tidytext)
library(dplyr)
library(stringr)
library(purrr)

## load twi

res = readRDS('/TRADE/all_SEX_TRADE_class_outliers.rds')
colnames(res) = c('twi','region','celltype')

## load perm twi significance

perm = readRDS('/TRADE/all_TRADE_perm_class.rds')

## load jk twi significance

p = readRDS('/TRADE/all_TRADE_jk_pval_class.rds')

## load twi loo range

loo = readRDS('~/TRADE/all_TRADE_loo_class.rds')

loo$sample = rownames(loo)
loo$sample = str_sub(loo$sample,1,5)
loopl = loo %>% group_by(region, celltype) %>% summarise(min = min(V1), max = max(V1))

## filter

meta = readRDS('meta_merged.rds')
check = reshape2::melt(meta[,c(3,12,41:151)]) 
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
res$key = paste(res$region, res$celltype)
res = subset(res, key %!in% remove$key)

## combine with perm p vals

perm$key = paste(perm$region, perm$celltype)
res = merge(res, perm[,c('key','p')], by = 'key')

## combine with loo p vals

p$key = paste(p$region, p$celltype)
res = merge(res, p[,c('key','p')], by = 'key')

## combine with loo range

loopl$key = paste(loopl$region, loopl$celltype)
res = merge(res, loopl[,c('key','min','max')], by = 'key')

## code for plotting

colnames(res)[c(5,6)] = c('p.perm','p.jk')

res$celltype = factor(res$celltype,levels = c(exn.levels, inn.levels, other.levels))
res$region = factor(res$region)
levels(res$region) = regions_short
res$key = paste(res$region, res$celltype)

res$perm.sig <- ifelse(res$p.perm < 0.05, 'TRUE', FALSE)

res$loo.sig <- ifelse(res$p.jk < 0.05, 'TRUE', FALSE)
res$loo.sig <- ifelse(res$twi > res$min & res$twi < res$max, res$loo.sig, FALSE)

table(res$p.perm < 0.05)
table(res$p.jk < 0.05)

ncol = which(c(exn.levels, inn.levels, other.levels) %!in% res$celltype)

## plot

res$max = ifelse(res$max > 0.4, 0.4, res$max)

# adjust p values

# res$p.perm.adj = p.adjust(res$p.perm, method = 'BH')
# res$p.jk.adj = p.adjust(res$p.jk, method = 'BH')
# res$perm.sig <- ifelse(res$p.perm.adj < 0.05, 'TRUE', FALSE)
# res$loo.sig <- ifelse(res$p.jk.adj < 0.05, 'TRUE', FALSE)
# res$loo.sig <- ifelse(res$twi > res$min & res$twi < res$max, res$loo.sig, FALSE)

# region x cell type 
# color = region

ggplot(res, aes(x = reorder(key, -twi), y = twi, fill = region)) +
  geom_bar(stat = 'identity', size = 1) +
  scale_fill_manual(values = region.colors) +
  theme_article() +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9)) +
  ylab("Transcriptome-wide impact") +
  geom_text(aes(label = ifelse(perm.sig, "*", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.8, size = 20 / .pt) +
  geom_text(aes(label = ifelse(loo.sig, "+", ""), group = key), 
            position = position_dodge(width = .9), vjust = -1.5, size = 10 / .pt) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = 'bottom',
        axis.title.x = element_blank(),
        legend.title = element_blank())

# within each cell type 
# color = region

ggplot(res, aes(x = reorder_within(x = region, by = -twi, within = celltype), y = twi, fill = region)) +
  geom_bar(stat = 'identity', size = 1) +
  scale_fill_manual(values = region.colors) +
  scale_x_reordered() +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label = ifelse(perm.sig, "*", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.4, size = 20 / .pt) +
  geom_text(aes(label = ifelse(loo.sig, "+", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.8, size = 10 / .pt) +
  facet_wrap(~celltype, scales = 'free_x') +
  ylab("Transcriptome-wide impact") +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = 'bottom',
        axis.title.x = element_blank(),
        legend.title = element_blank())

# region x cell type 
# color = cell type

ggplot(res, aes(x = reorder(key, -twi), y = twi, fill = celltype)) +
  geom_bar(stat = 'identity', size = 1) +
  scale_fill_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9),  alpha = 0.3) +
  geom_text(aes(label = ifelse(perm.sig, "*", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.4, size = 20 / .pt) +
  geom_text(aes(label = ifelse(loo.sig, "+", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.8, size = 10 / .pt) +
  ylab("Transcriptome-wide impact") +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = 'bottom',
        axis.title.x = element_blank(),
        legend.title = element_blank())

# within each region
# color = cell type

ggplot(res, aes(x = reorder_within(x = celltype, -twi, region), y = twi, fill = celltype)) +
  geom_bar(stat = 'identity', size = 1) +
  scale_fill_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) +
  ylab("Transcriptome-wide impact") +
  scale_x_reordered() +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9),  alpha = 0.3) +
  geom_text(aes(label = ifelse(perm.sig, "*", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.2, size = 20 / .pt) +
  geom_text(aes(label = ifelse(loo.sig, "+", ""), group = key), 
            position = position_dodge(width = .9), vjust = -0.6, size = 10 / .pt) +
  facet_wrap(~region, scales = 'free_x', ncol = 1) +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = 'none',
        axis.title.x = element_blank(),
        legend.title = element_blank())

# region x cell type 
# color = cell type
# outline = region

ggplot(res, aes(x = reorder(key, -twi), y = twi, fill = celltype, color = region)) +
  geom_bar(stat = 'identity', size = 1) +
  scale_color_manual(values = region.colors) +
  scale_fill_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(.9),  color = 'black',alpha = 0.3) +
  geom_text(aes(label = ifelse(perm.sig, "*", ""), group = key), 
            position = position_dodge(width = .9), color='black',vjust = -0.2, size = 20 / .pt) +
  geom_text(aes(label = ifelse(loo.sig, "+", ""), group = key), 
            position = position_dodge(width = .9), color='black',vjust = -0.6, size = 10 / .pt) +
  ylab("Transcriptome-wide impact") +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = 'none',
        axis.title.x = element_blank(),
        legend.title = element_blank())

## summarise

res %>% group_by(celltype) %>% summarise(mean = mean(twi, na.rm = T))
res %>% group_by(region) %>% summarise(mean = mean(twi, na.rm = T))

res$class = ifelse(res$celltype %in% exn.levels, 'exn', ifelse(res$celltype %in% inn.levels,'inn','other'))

mod = aov(twi ~ class, data = res)
summary(mod)
TukeyHSD(mod)

## plot summaries

ggplot(res, aes(x = reorder(celltype, -twi, FUN = median), y = twi, fill = celltype)) +
  geom_boxplot() +
  scale_fill_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title = element_blank(),
        legend.position = 'none')
  
ggplot(res, aes(x = reorder(region, -twi, FUN = median), y = twi, fill = region)) +
  geom_boxplot() +
  scale_fill_manual(values = region.colors) +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title = element_blank(),
        legend.position = 'none')

res2 = res
res2$class = ifelse(res2$celltype %in% exn.levels, 'exn', ifelse(res2$celltype %in% inn.levels,'inn','other'))
res2$class2 = ifelse(res2$celltype %in% other.levels, 'other', 'neuron')
levels(res2$region)
res2$region = factor(res2$region, levels = levels(res2$region)[c(6,1,4,3,2,5)])

ggplot(res2, aes(x = region, y = twi, fill = region)) +
  geom_boxplot() +
  scale_fill_manual(values = region.colors[c(6,1,4,3,2,5)]) +
  facet_wrap(~class, ncol = 1, scales = 'free_y') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title = element_blank(),
        legend.position = 'none')

mod = aov(twi ~ region, data = res2)
summary(mod)
TukeyHSD(mod)

mod = aov(twi ~ class + celltype + region, data = res2)
summary(mod)
TukeyHSD(mod)

mod = aov(twi ~ class + region, data = res2)
summary(mod)
TukeyHSD(mod)

mod = aov(twi ~ celltype + region, data = res2)
summary(mod)
TukeyHSD(mod)

mod = aov(twi ~ class2 + region, data = res2)
summary(mod)
TukeyHSD(mod)

mod = aov(twi ~ class2, data = res2)
summary(mod)
TukeyHSD(mod)

mod = aov(twi ~ region, data = subset(res2, class == 'other'))
mod = aov(twi ~ region, data = subset(res2, class == 'inn'))
mod = aov(twi ~ region, data = subset(res2, class == 'exn'))
summary(mod)
TukeyHSD(mod)

res2$region_code = ifelse(res2$region %in% c('AngGyr','RetCor'),'nb','biased')
mod = aov(twi ~ region_code, data = res2)
mod = aov(twi ~ region_code, data = subset(res2, class == 'other'))
mod = aov(twi ~ region_code, data = subset(res2, class == 'inn'))
mod = aov(twi ~ region_code, data = subset(res2, class == 'exn'))
summary(mod)

## compare TWI for sex and age

sextwi = readRDS('/TRADE/all_SEX_TRADE_class_outliers.rds')
colnames(sextwi) = c('twi','region','celltype')
sextwi$key = paste(sextwi$region, sextwi$celltype)
sextwi = subset(sextwi, key %!in% remove$key)
write.csv(sextwi, file = 'sextwi.csv')

agetwi = readRDS('/TRADE/all_AGE_TRADE_class_outliers.rds')
colnames(agetwi) = c('twi','region','celltype')
agetwi$key = paste(agetwi$region, agetwi$celltype)
agetwi = subset(agetwi, key %!in% remove$key)
write.csv(agetwi, file = 'agetwi.csv')

pl = data.frame(sextwi = sextwi$twi, agetwi = agetwi$twi, key = agetwi$key, region = agetwi$region, celltype = agetwi$celltype)
cor.test(pl$sextwi, pl$agetwi, method='spearman')
pl$celltype = factor(pl$celltype, levels = c(exn.levels, inn.levels, other.levels))

ggplot(pl, aes(x = sextwi, y = agetwi)) +
  geom_point(aes(x = sextwi, y = agetwi, color = celltype)) +
  #xlim(c(0,0.7)) + ylim(c(0,0.7)) +
  scale_color_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) +
  geom_smooth(method = 'lm') +
  theme_article()

ggplot(pl, aes(x = sextwi, y = agetwi)) +
  geom_point(aes(x = sextwi, y = agetwi, color = celltype)) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) +
  facet_wrap(~region, scales = 'free') +
  theme_article()

pl$class = ifelse(pl$celltype %in% exn.levels, 'exn', ifelse(pl$celltype %in% inn.levels,'inn','other'))

ggplot(pl, aes(x = sextwi, y = agetwi)) +
  geom_point(aes(x = sextwi, y = agetwi, color = celltype)) +
  scale_color_manual(values = c(exn.colors, inn.colors, other.colors)[-c(ncol)]) +
  facet_wrap(~class, scales = 'free') +
  geom_smooth(method = 'lm') +
  theme_article()

pl$diff = pl$agetwi - pl$sextwi

ggplot(pl, aes(x = celltype, y = diff)) +
  geom_boxplot() +
  theme_article()

## compare TWI to sex-biased GMV

gmv = readRDS('sexGMV.rds')

# gmv vs mean sexTWI across celltypes

res2 = res %>% group_by(region) %>% summarise(twi = mean(twi))

res_joined <- res2 %>%
  left_join(gmv, by = "region")

cor_results <- res_joined %>%
  summarise(
    cor_bias = cor(twi, bias, method = "spearman"),
    p_bias   = cor.test(twi, bias)$p.value,
    
    cor_absbias = cor(twi, absbias, method = "spearman"),
    p_absbias   = cor.test(twi, absbias)$p.value,
    
    cor_ranksigned = cor(twi, ranksigned, method = "spearman"),
    p_ranksigned   = cor.test(twi, ranksigned)$p.value,
    
    cor_rankabs = cor(twi, rankabs, method = "spearman"),
    p_rankabs   = cor.test(twi, rankabs)$p.value,
  )
View(cor_results)

ggplot(res_joined, aes(x = bias, y = twi)) +
  geom_point() + 
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_point(aes(x = bias, y = twi, color = region), size = 6) +
  geom_smooth(method='lm') +
  scale_color_manual(values = region.colors) +
  theme_article()

# gmv vs sexTWI within cell types

res_joined <- res %>%
  left_join(gmv, by = "region")

cor_results <- res_joined %>%
  group_by(celltype) %>%
  summarise(
    cor_bias = cor(twi, bias, method = "spearman"),
    p_bias   = cor.test(twi, bias, method = "spearman")$p.value,

    cor_absbias = cor(twi, absbias, method = "spearman"),
    p_absbias   = cor.test(twi, absbias, method = "spearman")$p.value,

    cor_ranksigned = cor(twi, ranksigned, method = "spearman"),
    p_ranksigned   = cor.test(twi, ranksigned, method = "spearman")$p.value,

    cor_rankabs = cor(twi, rankabs, method = "spearman"),
    p_rankabs   = cor.test(twi, rankabs, method = "spearman")$p.value,
  )
View(cor_results)
write.csv(cor_results, file = 'gmv-twi-correlations.csv')

cor_results$p_bias_adj = p.adjust(cor_results$p_bias)
cor_results$p_absbias_adj = p.adjust(cor_results$p_absbias)
cor_results$p_ranksigned_adj = p.adjust(cor_results$p_ranksigned)
cor_results$p_rankabs_adj = p.adjust(cor_results$p_rankabs)

library(tidyr)

cor_long <- cor_results %>%
  pivot_longer(
    cols = starts_with("cor_"),
    names_to = "metric",
    values_to = "correlation"
  )

pval_long <- cor_results %>%
  pivot_longer(
    cols = starts_with("p_"),    
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
    title = "Correlation of TWI with Bias Metrics by Cell Type",
    x = "Cell Type",
    y = "Pearson Correlation (r)",
    color = "Metric") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



