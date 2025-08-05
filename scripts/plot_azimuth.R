#!/usr/bin/env Rscript

# scripts/plot_azimuth.R
# Plot Azimuth annotation results across multiple reference mappings

library(ggplot2)
library(egg)
library(stringr)
library(dplyr)

##########
# choose one file to load
##########

c = readRDS("/azimuth/other_allen_M1.rds")
c = readRDS("/azimuth/other_allen_MCA.rds")
c = readRDS("/azimuth/other_BRAIN-BA40.rds")
c = readRDS("/azimuth/other_BRAIN-BA29-30.rds")
c = readRDS("/azimuth/other_BRAIN-BA5-7.rds")
c = readRDS("/azimuth/other_BRAIN-BA35-36.rds")
c = readRDS("/azimuth/other_BRAIN-gran-insula.rds")
c = readRDS("/azimuth/other_BRAIN-other.rds")

c = readRDS("/azimuth/inn_allen_M1.rds")
c = readRDS("/azimuth/inn_allen_MCA.rds")
c = readRDS("/azimuth/inn_BRAIN-BA40.rds")
c = readRDS("/azimuth/inn_BRAIN-BA29-30.rds")
c = readRDS("/azimuth/inn_BRAIN-BA5-7.rds")
c = readRDS("/azimuth/inn_BRAIN-BA35-36.rds")
c = readRDS("/azimuth/inn_BRAIN-gran-insula.rds")
c = readRDS("/azimuth/inn_BRAIN-inn.rds")

c = readRDS("/azimuth/exn_upper_allen_M1.rds")
c = readRDS("/azimuth/exn_upper_allen_MCA.rds")
c = readRDS("/azimuth/exn_upper_BRAIN-BA40.rds")
c = readRDS("/azimuth/exn_upper_BRAIN-BA29-30.rds")
c = readRDS("/azimuth/exn_upper_BRAIN-BA5-7.rds")
c = readRDS("/azimuth/exn_upper_BRAIN-BA35-36.rds")
c = readRDS("/azimuth/exn_upper_BRAIN-gran-insula.rds")
c = readRDS("/azimuth/exn_upper_BRAIN-exn.rds")

c = readRDS("/azimuth/exn_lower_allen_M1.rds")
c = readRDS("/azimuth/exn_lower_allen_MCA.rds")
c = readRDS("/azimuth/exn_lower_BRAIN-BA40.rds")
c = readRDS("/azimuth/exn_lower_BRAIN-BA29-30.rds")
c = readRDS("/azimuth/exn_lower_BRAIN-BA5-7.rds")
c = readRDS("/azimuth/exn_lower_BRAIN-BA35-36.rds")
c = readRDS("/azimuth/exn_lower_BRAIN-gran-insula.rds")
c = readRDS("/azimuth/exn_lower_BRAIN-exn.rds")


##########
# plot cell class
##########

propallen = prop.table(table(c$predicted.subclass_label, c$cell_class), 2) #allen
propallen = prop.table(table(c$predicted.supercluster_term, c$cell_class), 2) #brain region
propallen = prop.table(table(c$predicted.Subclass, c$cell_class), 2) #brain cell type
rownames(propallen)
colnames(propallen)

# filter other + allen M1
propallen = propallen[which(rownames(propallen) %in% c('Astro','Endo','Oligo','OPC','VLMC','Micro-PVM')),]
# filter other + allen MCA
propallen = propallen[which(rownames(propallen) %in% c('Astrocyte','Endothelial','Oligodendrocyte','OPC','Microglia','Pericyte','VLMC')),]
# filter other + brain region
propallen = propallen[which(rownames(propallen) %in% c('Astrocyte','Committed oligodendrocyte precursor','Fibroblast','Microglia','Oligodendrocyte','Vascular','Oligodendrocyte precursor')),]
# filter other + brain cell type
propallen = propallen[which(rownames(propallen) %in% c('Endothelial','VLMC','Microglia-PVM','Oligodendrocyte','OPC','Astrocyte')),]
# filter other
propallen = propallen[,which(colnames(propallen) != "Remove")]

# filter inn + allen M1
propallen = propallen[which(rownames(propallen) %in% c('Lamp5','Pvalb','Sncg','Sst','Vip')),]
# filter inn + allen MCA
propallen = propallen[which(rownames(propallen) %in% c('PVALB','SST','VIP','LAMP5','PAX6')),]
# filter inn + brain region
propallen = propallen[which(rownames(propallen) %in% c('CGE interneuron','MGE interneuron','LAMP5-LHX6 and Chandelier')),]
# filter inn + brain cell type
propallen = propallen[which(rownames(propallen) %in% c('Lamp5','Pvalb','Sncg','Sst','Vip','Chandelier','Lamp5 Lhx6','Pax6','Sst Chodl')),]

# filter exn + allen M1
propallen = propallen[which(rownames(propallen) %in% c('L2/3 IT','L5 IT')),]
propallen = propallen[which(rownames(propallen) %in% c('L5 IT','L5 ET','L5/6 NP','L6b','L6 CT','L6 IT','L6 IT Car3')),]
# filter exn + allen MCA
propallen = propallen[which(rownames(propallen) %in% c('L4 IT','L5 ET','L5/6 NP','IT')),]
propallen = propallen[which(rownames(propallen) %in% c('L5 ET','L6 IT','L6 IT Car3','L6b','L5/6 NP','IT','L6 CT')),]
# filter exn + brain region
propallen = propallen[which(rownames(propallen) %in% c('Upper-layer intratelencephalic','Deep-layer intratelencephalic','Deep-layer corticothalamic and 6b')),]
propallen = propallen[which(rownames(propallen) %in% c('Deep-layer near-projecting','Miscellaneous','Deep-layer intratelencephalic','Deep-layer corticothalamic and 6b')),]
# filter exn + brain cell type
propallen = propallen[which(rownames(propallen) %in% c('L2/3 IT','L4 IT')),]
propallen = propallen[which(rownames(propallen) %in% c('L5 IT','L5 ET','L5/6 NP','L6b','L6 CT','L6 IT','L6 IT Car3')),]

# plot

pheat_allen = slanter::sheatmap(propallen, color = viridis::viridis(100),
                  fontsize = 10,cellwidth = 10,cellheight = 10, cluster_rows = T, cluster_cols = F)

##########
# plot supertype
##########

propallen = prop.table(table(c$predicted.cluster_label, c$cell_class), 2) #allen
propallen = prop.table(table(c$predicted.supercluster_term, c$cell_class), 2) #brain region
propallen = prop.table(table(c$predicted.Supertype, c$cell_class), 2) #brain cell type

rownames(propallen)
colnames(propallen)

# filter other 
propallen = propallen[which(str_sub(rownames(propallen), 1, 3) %in% c('Ast','Oli','VLM','OPC','Per','End','Mic')),]
# filter inn 
propallen = propallen[which(str_sub(rownames(propallen), 1, 3) == 'Inh'),]
# filter inn 
propallen = propallen[which(str_sub(rownames(propallen), 1, 3) == 'Exc'),]
# filter all
propallen = propallen[which(rowSums(propallen) >= 0.2),]
propallen = propallen[which(rowSums(propallen) >= 0.1),] # inn BRAIN

pheat_allen = slanter::sheatmap(propallen, color = viridis::viridis(100),
                                fontsize = 6,cellwidth = 6,cellheight = 6, cluster_rows = T, cluster_cols = F)

