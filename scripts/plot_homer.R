library(ggplot2)
library(stringr)
library(biomaRt)

sexres = readRDS('sexres.rds')

hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
g2c = getBM(attributes = c('external_gene_name','chromosome_name'), mart = hsap)

##########
## homer
##########

homer = read.csv('homer-motifs.csv')
homer = unique(homer)
colnames(homer) = 'external_gene_name'
homer$external_gene_name <- toupper(homer$external_gene_name)
homer = merge(homer, g2c, by = 'external_gene_name', all.x = T)
View(subset(homer, chromosome_name %in% c('X','Y')))

homer = read.csv('/Users/decasienar/Desktop/Papers/snRNAseq/alevin_mapped_20k/homer/combined_homer_results.csv')

homer$p = 10^(homer$logp)

sh = c('ARE(NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer',
       'ERb(NR),IR3/Ovary-ERb-ChIP-Seq(GSE203391)/Homer',
       'ERE(NR),IR3/MCF7-ERa-ChIP-Seq(Unpublished)/Homer',
       'PGR(NR)/EndoStromal-PGR-ChIP-Seq(GSE69539)/Homer', # endometrium
       'PR(NR)/T47D-PR-ChIP-Seq(GSE31130)/Homer', # cell lines
       'ZFX(Zf)/mES-Zfx-ChIP-Seq(GSE11431)/Homer',
       'TFE3(bHLH)/MEF-TFE3-ChIP-Seq(GSE75757)/Homer',
       'ZBTB33(Zf)/GM12878-ZBTB33-ChIP-Seq(GSE32465)/Homer',
       'ZNF41(Zf)/HEK293-ZNF41.GFP-ChIP-Seq(GSE58341)/Homer',
       'ZNF711(Zf)/SHSY5Y-ZNF711-ChIP-Seq(GSE20673)/Homer',
       'Elf4(ETS)/BMDM-Elf4-ChIP-Seq(GSE88699)/Homer',
       'Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer',
       'Zic3(Zf)/mES-Zic3-ChIP-Seq(GSE37889)/Homer',
       'CDX4(Homeobox)/ZebrafishEmbryos-Cdx4.Myc-ChIP-Seq(GSE48254)/Homer',
       'Gata1(Zf)/K562-GATA1-ChIP-Seq(GSE18829)/Homer',
       'Sox3(HMG)/NPC-Sox3-ChIP-Seq(GSE33059)/Homer')

pl = subset(homer, motif %in% sh)
pl <- pl %>%
  group_by(cluster, input_type) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  ungroup()
pl <- pl %>%
  mutate(sig = case_when(
    padj < 0.05 & str_sub(input_type, -2) == "fn" ~ "fb",
    padj < 0.05 & str_sub(input_type, -2) == "mn" ~ "mb",
    padj < 0.05 & str_sub(input_type, -2) == "sn" ~ "sb",
    TRUE ~ "ns"
  ))
pl$cluster = as.factor(pl$cluster)
pl = subset(pl, input_type != 'homer-sn')
pl$plotp = ifelse(pl$input_type == 'homer-fn', -1*-log10(pl$padj), -log10(pl$padj))
table(pl$motif, pl$sig)
pl$motif = as.factor(pl$motif)
levels(pl$motif) = c('ARE','CDX4','ELF4','ELK1','ERB','ERA','GATA1','PGR1','PGR2','SOX3','TFE3','ZBTB33','ZFX','ZIC3','ZNF41','ZNF711')
unique(pl$motif[which(pl$motif %!in% sexres$gene)]) # CDX4  SOX3  PGR1  ERE   ERB   PGR2  GATA1
pl = subset(pl, motif %!in% c('CDX4',  'SOX3', 'GATA1'))
pl$motif = droplevels(pl$motif)
pl$motif = factor(pl$motif, levels = c('ARE','ERB','ERA','PGR1','PGR2','ELF4','ELK1','TFE3','ZBTB33','ZFX','ZIC3','ZNF41','ZNF711'))
pl$motif = factor(pl$motif, levels = rev(levels(pl$motif)))
pl$code = ifelse(pl$p < 0.05, "+", NA)
pl$code = ifelse(pl$padj < 0.1, "*", pl$code)
range(pl$plotp)
pl$target_perc <- as.numeric(sub("%", "", pl$target_perc)) / 100
pl$back_perc <- as.numeric(sub("%", "", pl$back_perc)) / 100
pl$or = pl$target_perc / pl$back_perc
pl$plotor = ifelse(pl$input_type == 'homer-fn', -1*pl$or, pl$or)
range(pl$plotor)
pl$group = ifelse(pl$input_type == 'homer-fn', 'female-biased', 'male-biased')
#pl$s = ifelse(pl$motif %in% c('ARE','ERB','ERA','PGR1','PGR2'),'hormone','Xchr')
pl$s = ifelse(pl$motif %in% c('ARE'),'AR','Xchr')
pl$s = ifelse(pl$motif %in% c('ERA','ERB'),'ER',pl$s)
pl$s = ifelse(pl$motif %in% c('PGR1','PGR2'),'PR',pl$s)
table(pl$s)

ggplot(pl, aes(x = cluster, y = motif, fill = plotp)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-5,5), low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(s~group, scales = 'free', space = 'free', drop = T) +
  theme_article() +
  geom_text(label = pl$code, size = 5) +
  theme(axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 18),
        axis.title = element_blank())

###########
## HOCOMOCO
###########

HOCOMOCO = read.csv('HOCOMOCO.csv')
colnames(HOCOMOCO) = 'external_gene_name'
HOCOMOCO = unique(HOCOMOCO)
HOCOMOCO = merge(HOCOMOCO, g2c, by = 'external_gene_name', all.x = T)
View(subset(HOCOMOCO, chromosome_name %in% c('X','Y')))
View(subset(HOCOMOCO, external_gene_name == 'SRY'))

HOCOMOCO = read.csv('/Users/decasienar/Desktop/Papers/snRNAseq/alevin_mapped_20k/homer/combined_homer_HOCOMOCO_results.csv')
HOCOMOCO$p = 10^(HOCOMOCO$logp)

sh = c('ANDR_HUMAN.H11MO.0.A',
       'ESR1_HUMAN.H11MO.0.A',
       'ESR2_HUMAN.H11MO.0.A',
       'PRGR_HUMAN.H11MO.0.A', 
       'FOXP3_HUMAN.H11MO.0.D',
       'GATA1_HUMAN.H11MO.1.A',
       'HSFY1_HUMAN.H11MO.0.D',
       'SOX3_HUMAN.H11MO.0.B',
       'SHOX_HUMAN.H11MO.0.D',
       'ZFX_HUMAN.H11MO.0.A',
       'SRY_HUMAN.H11MO.0.B',
       'ARX_HUMAN.H11MO.0.D',
       'ELK1_HUMAN.H11MO.0.B',
       'FOXO4_HUMAN.H11MO.0.C',
       'KLF8_HUMAN.H11MO.0.C',
       'MECP2_HUMAN.H11MO.0.C',
       'NR0B1_HUMAN.H11MO.0.D',
       'TAF1_HUMAN.H11MO.0.A',
       'TFE3_HUMAN.H11MO.0.B',
       'ZIC3_HUMAN.H11MO.0.B',
       'ZNF41_HUMAN.H11MO.0.C',
       'ZBED1_HUMAN.H11MO.0.D')

sh[which(sh %!in% HOCOMOCO$motif)]

pl2 = subset(HOCOMOCO, motif %in% sh)
pl2 <- pl2 %>%
  group_by(cluster, input_type) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  ungroup()
pl2 <- pl2 %>%
  mutate(sig = case_when(
    padj < 0.05 & str_sub(input_type, -11) == "fn-HOCOMOCO" ~ "fb",
    padj < 0.05 & str_sub(input_type, -11) == "mn-HOCOMOCO" ~ "mb",
    padj < 0.05 & str_sub(input_type, -11) == "sn-HOCOMOCO" ~ "sb",
    TRUE ~ "ns"
  ))
pl2$cluster = as.factor(pl2$cluster)
pl2 = subset(pl2, input_type != 'homer-sn-HOCOMOCO')
pl2$plotp = ifelse(pl2$input_type == 'homer-fn-HOCOMOCO', -1*-log10(pl2$padj), -log10(pl2$padj))
table(pl2$motif, pl2$sig)
pl2$motif = as.factor(pl2$motif)
levels(pl2$motif) = c('AR','ARX','ELK1','ERA','ERB','FOXO4','FOXP3','GATA1','HSFY1','KLF8','MECP2',
                    'NR0B1','PGR','SHOX','SOX3','SRY','TAF1','TFE3','ZBED1','ZFX','ZIC3','ZNF41')
pl2$motif = factor(pl2$motif, levels = c('AR','ERA','ERB','PGR','ARX','ELK1','FOXO4','FOXP3','GATA1','HSFY1','KLF8','MECP2',
                                       'NR0B1','SOX3','TAF1','TFE3','ZFX','ZIC3','ZNF41','SHOX','ZBED1','SRY'))
unique(pl2$motif[which(pl2$motif %!in% sexres$gene)]) # SHOX  GATA1 SOX3  HSFY1 ERA   FOXP3 ERB   SRY 
pl2 = subset(pl2, motif %!in% c('SHOX','GATA1','SOX3','HSFY1','FOXP3','SRY'))
pl2$motif = droplevels(pl2$motif)
pl2$motif = factor(pl2$motif, levels = rev(levels(pl2$motif)))
pl2$code = ifelse(pl2$p < 0.05, "+", NA)
pl2$code = ifelse(pl2$padj < 0.1, "*", pl2$code)
range(pl2$plotp)
pl2$target_perc <- as.numeric(sub("%", "", pl2$target_perc)) / 100
pl2$back_perc <- as.numeric(sub("%", "", pl2$back_perc)) / 100
pl2$or = pl2$target_perc / pl2$back_perc
pl2$plotor = ifelse(pl2$input_type == 'homer-fn-HOCOMOCO', -1*pl2$or, pl2$or)
range(pl2$plotor)
pl2$group = ifelse(pl2$input_type == 'homer-fn-HOCOMOCO', 'female-biased', 'male-biased')
#pl2$s = ifelse(pl2$motif %in% c('AR','ERA','ERB','PGR'),'hormone','Xchr')
#pl2$s = ifelse(pl2$motif %in% c('ZBED1'),'PAR',pl2$s)
pl2$s = ifelse(pl2$motif %in% c('AR'),'AR','Xchr')
pl2$s = ifelse(pl2$motif %in% c('ERA','ERB'),'ER',pl2$s)
pl2$s = ifelse(pl2$motif %in% c('PGR'),'PR',pl2$s)
pl2$s = ifelse(pl2$motif %in% c('ZBED1'),'PAR',pl2$s)
table(pl2$s)

ggplot(pl2, aes(x = cluster, y = motif, fill = plotp)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-5,5), low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(s~group, scales = 'free', space = 'free') +
  theme_article() +
  geom_text(label = pl2$code, size = 5) +
  theme(axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 18),
        axis.title = element_blank())

##############
## combined
##############

c = rbind(cbind(pl, source = 'homer'), cbind(pl2, source = 'hocomoco'))
c$mot = paste(c$motif, c$source)
c$mot = as.factor(c$mot)
levels(c$mot)
c$mot = factor(c$mot, levels = rev(levels(c$mot)))
c$s = as.factor(c$s)
c$s = factor(c$s, levels = c('AR','ER','PR','PAR','Xchr'))
write.csv(c, 'tf-results.csv')

ggplot(c, aes(x = cluster, y = mot, fill = plotp)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-5,5), low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = 'white', na.value = "grey90") +
  facet_grid(s~group, scales = 'free', space = 'free') +
  theme_article() +
  geom_text(label = c$code, size = 5, vjust = 'center') +
  theme(axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 18),
        axis.title = element_blank())
