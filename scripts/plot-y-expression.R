library(dplyr)
library(ggplot2)
library(egg)
library(stringr)

pseudo = readRDS('/pseudobulk/all_res_normalized_expression_class_corrected_outliers.rds')
gene2chrom = readRDS('gene2chrom_mash.rds')
meta = readRDS('meta_merged.rds')

# XY males only
# Y chr genes only

m = subset(meta, sex == 'Male')
pseudo = pseudo[which(rownames(pseudo) %in% subset(gene2chrom, chromosome_name == 'Y')$external_gene_name),which(str_sub(colnames(pseudo),1,5) %in% m$new.id)]

pseudo$gene = rownames(pseudo)
pl = reshape2::melt(pseudo)
d = str_split_fixed(pl$variable, pattern = '\\.', 3)
pl = cbind(pl, d)
colnames(pl) = c('gene','key','log2cpm','sample','region','celltype')
pp = pl %>% group_by(gene, region, celltype) %>% summarise(mean = mean(log2cpm, na.rm = T))
pp$cpm = 2^pp$mean
pp$celltype = factor(pp$celltype, levels = c(exn.levels, inn.levels, other.levels))

# add missing gametologs

gams = c('PRKX','PRKY','ZFX','ZFY','USP9X','USP9Y','TBL1X','TBL1Y','NLGN4Y','NLGN4X','UTY','KDM6A','TMSB4X','TMSB4Y','TGIF2LX','TGIF2LY',
         'RPS4X','RPS4Y1','DDX3Y','DDX3X','EIF1AY','EIF1AX','KDM5C','KDM5D','PCDH11Y','PCDH11X','TXLNG','TXLNGY','AMELX','AMELY','SOX3','SRY')
pp$gam = ifelse(pp$gene %in% gams, 'gam', 'x')
gams[which(gams %!in% pp$gene)]
pp = rbind(pp, data.frame(gene = 'TMSB4Y', region = 'Angular_Gyrus', celltype = 'Oligo', mean = NA, cpm = NA, gam = 'gam'))
pp = rbind(pp, data.frame(gene = 'TGIF2LY', region = 'Angular_Gyrus', celltype = 'Oligo', mean = NA, cpm = NA, gam = 'gam'))
pp = rbind(pp, data.frame(gene = 'TXLNGY', region = 'Angular_Gyrus', celltype = 'Oligo', mean = NA, cpm = NA, gam = 'gam'))
pp = rbind(pp, data.frame(gene = 'AMELY', region = 'Angular_Gyrus', celltype = 'Oligo', mean = NA, cpm = NA, gam = 'gam'))
pp = rbind(pp, data.frame(gene = 'SRY', region = 'Angular_Gyrus', celltype = 'Oligo', mean = NA, cpm = NA, gam = 'gam'))
pp = tibble(pp)
pp = pp %>% complete(gene, region, celltype)
pp$gam = ifelse(pp$gene %in% gams, 'gam', 'x')
pp$celltype = factor(pp$celltype, levels = c(exn.levels, inn.levels, other.levels))

## filter

check = reshape2::melt(meta[,c(3,4,44:154)]) 
check = data.frame(table(check$region, check$variable, check$value > 50))
check = subset(check, Var3 == TRUE)
check$key = paste(check$Var1, check$Var2)
remove = subset(check, Freq < 10)
pp$key = paste(pp$region, pp$celltype)
pp = subset(pp, key %!in% remove$key)
pp$region = factor(pp$region)
levels(pp$region) = regions_short

## plot

ggplot(pp, aes(x = celltype, y = gene, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey50") +
  facet_grid(gam~region, scales = 'free', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        strip.text = element_text(size = 18),
        axis.title = element_blank())

## summarise

p = pp %>% group_by(gene, celltype) %>% summarise(sum = sum(mean, na.rm = T))
View(p)

p = pp %>% group_by(region, celltype) %>% summarise(sum = sum(mean, na.rm = T))
View(p)

p = pp[complete.cases(pp),] %>% group_by(region, celltype) %>% summarise(n = n())
View(p)
