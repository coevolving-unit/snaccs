library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(ggplot2)
library(egg)
library(pheatmap)
library(tidytext)
library(stringr)
library(forcats)
library(reshape2)

#### load data ####

tx2gene = readRDS('gene2chrom_mash.rds')
sig = read.csv('significant_genes_per_cluster.csv')
bg = read.csv('background_genes_per_cluster.csv') 

# autosomes only
sig = subset(sig, gene %in% subset(tx2gene, chromosome_name %in% c(1:22))$external_gene_name)
bg = subset(bg, gene %in% subset(tx2gene, chromosome_name %in% c(1:22))$external_gene_name)

table(sig$bias, sig$cluster)

#### end ####

#### run GO ####

clus = unique(sig$cluster)
bp_go = data.frame()
cc_go = data.frame()

for(i in 1:length(clus)){
  
  print(clus[i])
  mn = subset(sig, bias == 'male' & cluster == clus[i])$gene
  fn = subset(sig, bias == 'female' & cluster == clus[i])$gene
  bn = subset(bg, cluster == clus[i])$gene
  
  m_go <- enrichGO(gene = mn,universe = bn,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 1)
  m_go = m_go@result
  rownames(m_go) = NULL
  m_go$cluster = i
  m_go$bias = 'male'
  bp_go = rbind(bp_go, m_go)

  f_go <- enrichGO(gene = fn,universe = bn,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 1)
  f_go = f_go@result
  rownames(f_go) = NULL
  f_go$cluster = i
  f_go$bias = 'female'
  bp_go = rbind(bp_go, f_go)

  m_go <- enrichGO(gene = mn,universe = bn,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "CC",pAdjustMethod = "BH", pvalueCutoff = 1)
  m_go = m_go@result
  rownames(m_go) = NULL
  m_go$cluster = i
  m_go$bias = 'male'
  cc_go = rbind(cc_go, m_go)

  f_go <- enrichGO(gene = fn,universe = bn,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "CC",pAdjustMethod = "BH", pvalueCutoff = 1)
  f_go = f_go@result
  rownames(f_go) = NULL
  f_go$cluster = i
  f_go$bias = 'female'
  cc_go = rbind(cc_go, f_go)
  
}

saveRDS(bp_go, file = 'bp_go_auto.rds')
saveRDS(cc_go, file = 'cc_go_auto.rds')

#### end ####

#### plot ####

## biological processes

bp_go = readRDS('bp_go_auto.rds')
bp_go = subset(bp_go, pvalue < 0.05)

bp_go = bp_go[order(bp_go$FoldEnrichment, decreasing = F),]
pn = bp_go %>% group_by(cluster, bias) %>% dplyr::slice(1:3) 
pn$FoldEnrichment = ifelse(pn$bias == 'male', pn$FoldEnrichment, -1*pn$FoldEnrichment)
ggplot(pn, aes(y = reorder_within(Description, FoldEnrichment, cluster), x = FoldEnrichment, fill = bias)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = sex.colors) +
  facet_grid(cluster~., scales = 'free') +
  theme_article() +
  scale_y_reordered() + 
  theme(axis.title.y = element_blank(),
        legend.position = 'none')

pm = subset(bp_go, Description %in% pn$Description)
pm$FoldEnrichment = ifelse(pm$bias == 'male', pm$FoldEnrichment, -1*pm$FoldEnrichment)
pm$cluster = as.factor(pm$cluster)
pm$Description = as.factor(pm$Description)
View(pm %>% group_by(Description, cluster) %>% summarise(n=n()))
cb = data.frame(pm %>% group_by(Description, cluster) %>% summarise(n=n()))
table(cb$n)

pm = merge(cb, pm[,c('Description','cluster','bias')], by = c('Description','cluster'))
pm$FoldEnrichment = ifelse(pm$n == 2, '0', ifelse(pm$n == 1 & pm$bias == 'female', -1, 1))
pm$FoldEnrichment = as.numeric(pm$FoldEnrichment)
table(pm$FoldEnrichment)

length(unique(pm$Description))

unique(pm$Description)

functional_groups <- list(
  "Development/Morphogenesis" = c(
    "axon development", "axonogenesis", "cell growth", 
    "cellular anatomical entity morphogenesis", "developmental process involved in reproduction",
    "multicellular organismal reproductive process", "nuclear division", 
    "positive regulation of cell population proliferation","negative regulation of cell population proliferation", "regulation of cell projection organization",
    "vasculature development", "blood vessel development", "sensory system development",
    "gamete generation", "spermatogenesis", "sexual reproduction", "wound healing"
  ),
  
  "Signaling/Communication" = c(
    "cell-cell signaling by wnt", "G protein-coupled receptor signaling pathway", 
    "regulation of canonical Wnt signaling pathway", "regulation of trans-synaptic signaling",
    "regulation of establishment of protein localization", "regulation of hydrolase activity",
    "response to growth factor", "response to inorganic substance", 
    "trans-synaptic signaling", "Wnt signaling pathway", 
    "transmembrane receptor protein tyrosine kinase signaling pathway", 
    "small GTPase-mediated signal transduction", "positive regulation of catalytic activity"
  ),
  
  "Neural/Synaptic" = c(
    "chemical synaptic transmission", "synaptic signaling", 
    "synapse organization", "secretion by cell"
  ),
  
  "Cytoskeleton/Motility" = c(
    "actin cytoskeleton organization", "chromosome organization", 
    "cell junction organization", "mitochondrion organization", 
    "regulation of cytoskeleton organization", 
    "negative regulation of cellular component organization",
    "positive regulation of cell motility", "positive regulation of cell migration"
  ),
  
  "Metabolism/Transport" = c(
    "lipid biosynthetic process", "monocarboxylic acid metabolic process",
    "organophosphate biosynthetic process", "protein maturation", 
    "carbohydrate derivative biosynthetic process", 
    "cellular modified amino acid metabolic process",
    "inorganic cation transmembrane transport", 
    "inorganic ion transmembrane transport",
    "monoatomic cation transmembrane transport", 
    "multicellular organismal-level homeostasis"
  ),
  
  "Apoptosis/Immune" = c(
    "apoptotic signaling pathway", "intrinsic apoptotic signaling pathway", 
    "DNA repair", "defense response to symbiont", "immune response-regulating cell surface receptor signaling pathway", 
    "viral life cycle", "viral process"
  )
)

table(table(reshape2::melt(functional_groups)$value))
unique(pm$Description)[which(unique(pm$Description) %!in% melt(functional_groups)$value)]
pm$FunctionalGroup <- sapply(pm$Description, function(desc) {
  group_name <- names(which(sapply(functional_groups, function(terms) desc %in% terms)))
  if (length(group_name) == 0) return(NA) else return(group_name)})
table(pm$FunctionalGroup)

pm$Description = gsub("development", "dev", pm$Description)
pm$Description = gsub("transmembrane receptor ", "", pm$Description)
pm$Description = gsub("positive", "pos", pm$Description)
pm$Description = gsub("negative", "neg", pm$Description)
pm$Description = gsub("regulation", "reg", pm$Description)
pm$Description = gsub("regulating", "reg", pm$Description)
pm$Description = gsub("signaling", "sig", pm$Description)
pm$Description = gsub("response", "resp", pm$Description)

unique(pm$Description)

ggplot(pm, aes(x = cluster, y = Description, fill = FoldEnrichment)) +
  geom_tile() +
  facet_grid(FunctionalGroup~., space = 'free', scales = 'free', switch = "y") +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = '#8E89A0', na.value = "grey90") +
  theme_article() +
  theme(strip.placement = "outside", axis.title = element_blank())

pm$FunctionalGroupColor <- recode(pm$FunctionalGroup, 
                                  "Neuronal Development & Synaptic Signaling" = "#F0E422",  
                                  "Development & Morphogenesis" = "#FF7F00",         
                                  "Cell Adhesion, Cytoskeleton, & Motility" = "#4DAF4A",
                                  "Immune Response & Defense" = "#377EB8",  
                                  "Reproduction & Developmental Processes" = "#E41A1C", 
                                  "Transport, Biosynthesis & Metabolic Regulation" = "#9B59B6")

color_code = unique(pm[,c('Description','FunctionalGroupColor')])
color_code = color_code[match(levels(pm$Description), color_code$Description),]
ggplot(pm, aes(x = cluster, y = Description, fill = FoldEnrichment)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = '#8E89A0', na.value = "grey90") +
  geom_point(aes(x = cluster, y = Description, color = FunctionalGroup), size = 4, alpha = 0) +
  scale_color_manual(values = setNames(pm$FunctionalGroupColor, pm$FunctionalGroup)) +
  theme_article() +
  theme(axis.text.y = element_text(color = color_code$FunctionalGroupColor, size = 10),  
    axis.ticks.y = element_line(color = "black"), 
    axis.line.y = element_line(color = "black")) +
  guides(fill = guide_colorbar(title = "Fold Enrichment"), 
    color = guide_legend(title = "Functional Group", override.aes = list(size = 4, shape = 16, alpha = 1)))

# cellular compartments
                                   
cc_go = readRDS('cc_go_auto.rds')
cc_go = subset(cc_go, pvalue < 0.05)

cc_go = cc_go[order(cc_go$FoldEnrichment, decreasing = F),]
pn = cc_go %>% group_by(cluster, bias) %>% dplyr::slice(1:3) 
pn$FoldEnrichment = ifelse(pn$bias == 'male', pn$FoldEnrichment, -1*pn$FoldEnrichment)
ggplot(pn, aes(y = reorder_within(Description, FoldEnrichment, cluster), x = FoldEnrichment, fill = bias)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = sex.colors) +
  facet_grid(cluster~., scales = 'free') +
  theme_article() +
  scale_y_reordered() + 
  theme(axis.title.y = element_blank(),
        legend.position = 'none')

pm = subset(cc_go, Description %in% pn$Description)
pm$FoldEnrichment = ifelse(pm$bias == 'male', pm$FoldEnrichment, -1*pm$FoldEnrichment)
pm$cluster = as.factor(pm$cluster)
pm$Description = as.factor(pm$Description)
View(pm %>% group_by(Description, cluster) %>% summarise(n=n()))
cb = data.frame(pm %>% group_by(Description, cluster) %>% summarise(n=n()))
table(cb$n)

# when F and M enriched, code both
pm = merge(cb, pm[,c('Description','cluster','bias')], by = c('Description','cluster'))
pm$FoldEnrichment = ifelse(pm$n == 2, '0', ifelse(pm$n == 1 & pm$bias == 'female', -1, 1))
pm$FoldEnrichment = as.numeric(pm$FoldEnrichment)
table(pm$FoldEnrichment)

unique(pm$Description)

functional_groups <- list(
  "Synapse" = c(
    "axon", 
    "glutamatergic synapse", 
    "neuron to neuron synapse", 
    "postsynaptic density", 
    "postsynaptic specialization", 
    "presynapse"
  ),
  
  "Plasma Membrane" = c(
    "adherens junction", 
    "apical part of cell", 
    "apical plasma membrane", 
    "basal part of cell", 
    "basolateral plasma membrane", 
    "cell projection membrane", 
    "cell surface", 
    "cell-substrate junction", 
    "cytoplasmic side of plasma membrane", 
    "plasma membrane protein complex", 
    "side of membrane", 
    "receptor complex"
  ),
  
  "Mitochondria" = c(
    "mitochondrial envelope", 
    "mitochondrial inner membrane", 
    "mitochondrial intermembrane space", 
    "mitochondrial membrane", 
    "mitochondrial outer membrane", 
    "mitochondrial protein-containing complex", 
    "organelle inner membrane"
  ),
  
  "Vesicle Trafficking" = c(
    "endoplasmic reticulum lumen", 
    "endoplasmic reticulum-Golgi intermediate compartment", 
    "Golgi membrane", 
    "late endosome", 
    "late endosome membrane", 
    "recycling endosome", 
    "transport vesicle"
  ),
  
  "Extracellular" = c(
    "collagen-containing extracellular matrix", 
    "external encapsulating structure", 
    "extracellular matrix"
  ),
  
  "Nuclear Structures" = c(
    "centriole", 
    "condensed nuclear chromosome", 
    "methyltransferase complex", 
    "supramolecular fiber", 
    "supramolecular polymer", 
    "transferase complex, transferring phosphorus-containing groups", 
    "transmembrane transporter complex", 
    "transporter complex"
  )
)

unique(pm$Description)[which(unique(pm$Description) %!in% melt(functional_groups)$value)]
pm$FunctionalGroup <- sapply(pm$Description, function(desc) {
  group_name <- names(which(sapply(functional_groups, function(terms) desc %in% terms)))
  if (length(group_name) == 0) return(NA) else return(group_name)})
table(pm$FunctionalGroup)

unique(pm$Description)
pm$Description = gsub(", transferring phosphorus-containing groups", "", pm$Description)
pm$Description = gsub("endoplasmic reticulum", "ER", pm$Description)

ggplot(pm, aes(x = cluster, y = Description, fill = FoldEnrichment)) +
  geom_tile() +
  facet_grid(FunctionalGroup~., space = 'free', scales = 'free', switch = "y") +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = '#8E89A0', na.value = "grey90") +
  theme_article() +
  theme(strip.placement = "outside", axis.title = element_blank())

pm$FunctionalGroupColor <- recode(pm$FunctionalGroup, 
                                  "Membranes and Extracellular Components" = "#1F77B4",  # Blue
                                  "Cytoskeleton and Related Structures" = "#FF7F0E",      # Orange
                                  "Synaptic and Neuronal Compartments" = "#2CA02C",       # Green
                                  "Organelles and Endosomes" = "#D62728",                 # Red
                                  "Cell Junctions and Membrane Bound Complexes" = "#9467BD")  # Purple

color_code = unique(pm[,c('Description','FunctionalGroupColor')])
color_code = color_code[match(levels(pm$Description), color_code$Description),]
ggplot(pm, aes(x = cluster, y = Description, fill = FoldEnrichment)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[1], high = sex.colors[2], midpoint = 0, mid = '#8E89A0', na.value = "grey90") +
  geom_point(aes(x = cluster, y = Description, color = FunctionalGroup), size = 4, alpha = 0) +
  scale_color_manual(values = setNames(pm$FunctionalGroupColor, pm$FunctionalGroup)) +
  theme_article() +
  theme(axis.text.y = element_text(color = color_code$FunctionalGroupColor, size = 10),  
        axis.ticks.y = element_line(color = "black"), 
        axis.line.y = element_line(color = "black")) +
  guides(fill = guide_colorbar(title = "Fold Enrichment"), 
         color = guide_legend(title = "Functional Group", override.aes = list(size = 4, shape = 16, alpha = 1)))

#### end ####
