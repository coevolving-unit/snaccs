library(tidyr)
library(dplyr)

`%!in%` = Negate(`%in%`)

## load data

gene = read.csv('genotype/genedat.csv')
gene = gene[,-1]
cell = read.csv('genotype/celldat.csv')
cell = cell[,-1]

gene = subset(gene, tot_n_snps > 2)

combined <- gene[,c("orig.ident", "Domain.category", "Gene", "Annotation", "tot_n_snps", "bi_n_snps", "sex")] %>%
  left_join(
    cell[,c("orig.ident", "Gene", "Annotation", "total_cells", "cells_with_bi", "sex")],
    by = c("orig.ident", "Gene", "Annotation", "sex")
  )
combined$Domain.category = ifelse(combined$Domain.category == 'PAR', 'PAR', 'nonPAR')
combined$Domain.category = ifelse(combined$Gene == 'LINC00106', 'PAR', combined$Domain.category)
combined$Domain.category = ifelse(combined$Gene == 'LINC00685', 'PAR', combined$Domain.category)
combined$Domain.category[is.na(combined$Domain.category)] = 'nonPAR'

h = combined %>% group_by(Gene, sex, Domain.category) %>% summarise(sum = sum(cells_with_bi))
gene_summary <- h %>%
  group_by(Gene, Domain.category) %>%
  summarize(female_sum = sum(sum[sex == "Female"], na.rm = TRUE),
    male_sum = sum(sum[sex == "Male"], na.rm = TRUE))
gene_summary = subset(gene_summary, Domain.category != 'PAR')
#gene_summary = subset(gene_summary, Domain.category == 'PAR')
gene_summary = subset(gene_summary, !(female_sum == 0 & male_sum == 0))
gene_summary = subset(gene_summary, female_sum != 0)
nrow(gene_summary) 
both_female_and_male_positive <- gene_summary %>%
  filter(female_sum > 0 & male_sum > 0) %>%
  nrow()
both_female_and_male_positive 
female_positive_male_zero <- gene_summary %>%
  filter(female_sum > 0 & male_sum == 0) %>%
  nrow()
female_positive_male_zero 

h = combined %>% group_by(Gene, Annotation, sex, Domain.category) %>% summarise(sum = sum(cells_with_bi))
gene_summary <- h %>%
  group_by(Annotation, Gene, Domain.category) %>%
  summarize(female_sum = sum(sum[sex == "Female"], na.rm = TRUE),
            male_sum = sum(sum[sex == "Male"], na.rm = TRUE))
gene_summary = subset(gene_summary, Domain.category != 'PAR')
#gene_summary = subset(gene_summary, Domain.category == 'PAR')
gene_summary = subset(gene_summary, !(female_sum == 0 & male_sum == 0))
gene_summary = subset(gene_summary, female_sum != 0)
nrow(gene_summary) 
both_female_and_male_positive <- gene_summary %>%
  filter(female_sum > 0 & male_sum > 0) %>%
  nrow()
both_female_and_male_positive  
female_positive_male_zero <- gene_summary %>%
  filter(female_sum > 0 & male_sum == 0) %>%
  nrow()
female_positive_male_zero 
female_positive_male_1_0 <- gene_summary %>%
  filter(female_sum > 0 & male_sum <= 1) %>%
  nrow()
female_positive_male_1_0 

remove = c('FRMPD4','DMD','PTCHD1-AS','IDS','TBL1X','DANT2','FTX','IL1RAPL1','LINC00632')
combined = subset(combined, Gene %!in% remove)
combined$orig.ident = as.character(combined$orig.ident)

summary_df <- combined %>%
  group_by(Gene, Domain.category, sex) %>%
  summarise(
    sum_tot_n_snps = sum(tot_n_snps, na.rm = TRUE),
    sum_bi_n_snps = sum(bi_n_snps, na.rm = TRUE),
    sum_total_cells = sum(total_cells, na.rm = TRUE),
    sum_cells_with_bi = sum(cells_with_bi, na.rm = TRUE),
    count_individuals = length(unique(orig.ident)),
    .groups = 'drop'
  )

summary_wide <- summary_df %>%
  pivot_wider(
    names_from = sex,
    values_from = c(sum_tot_n_snps, sum_bi_n_snps, sum_total_cells, sum_cells_with_bi, count_individuals),
    names_sep = "."
  ) %>%
  dplyr::select(
    Gene, Domain.category,
    starts_with("sum_tot_n_snps.Female"), starts_with("sum_bi_n_snps.Female"),
    starts_with("sum_total_cells.Female"), starts_with("sum_cells_with_bi.Female"),
    starts_with("count_individuals.Female"),
    starts_with("sum_tot_n_snps.Male"), starts_with("sum_bi_n_snps.Male"),
    starts_with("sum_total_cells.Male"), starts_with("sum_cells_with_bi.Male"),
    starts_with("count_individuals.Male")
  )

summary_wide <- subset(summary_wide, !(sum_bi_n_snps.Female == 0 & sum_bi_n_snps.Male == 0))

# summary_wide$fbr = summary_wide$sum_bi_n_snps.Female / summary_wide$sum_tot_n_snps.Female
# summary_wide$mbr = summary_wide$sum_bi_n_snps.Male / summary_wide$sum_tot_n_snps.Male
# summary_wide$check = ifelse(summary_wide$fbr > summary_wide$mbr, 'ok', 'check')
summary_wide$check = ifelse(summary_wide$sum_cells_with_bi.Female > summary_wide$sum_cells_with_bi.Male, 'ok', 'check')

summary_wide = subset(summary_wide, !(check == 'check' & Domain.category != 'PAR'))
View(summary_wide)

summary_df <- combined %>%
  group_by(Gene, Domain.category, Annotation, sex) %>%
  summarise(
    sum_tot_n_snps = sum(tot_n_snps, na.rm = TRUE),
    sum_bi_n_snps = sum(bi_n_snps, na.rm = TRUE),
    sum_total_cells = sum(total_cells, na.rm = TRUE),
    sum_cells_with_bi = sum(cells_with_bi, na.rm = TRUE),
    count_bi_individuals = n_distinct(orig.ident[cells_with_bi > 0]),
    count_individuals = length(unique(orig.ident)),
    .groups = 'drop'
  )

summary_wide <- summary_df %>%
  group_by(Gene, Domain.category, Annotation, sex) %>%
  summarise(
    sum_tot_n_snps = sum(sum_tot_n_snps, na.rm = TRUE),
    sum_bi_n_snps = sum(sum_bi_n_snps, na.rm = TRUE),
    sum_total_cells = sum(sum_total_cells, na.rm = TRUE),
    sum_cells_with_bi = sum(sum_cells_with_bi, na.rm = TRUE),
    count_bi_individuals = sum(count_bi_individuals, na.rm = TRUE),
    count_individuals = sum(count_individuals, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(sum_tot_n_snps, sum_bi_n_snps, sum_total_cells, sum_cells_with_bi, count_bi_individuals,count_individuals),
    names_sep = "."
  ) %>%
  dplyr::select(
    Gene, Domain.category, Annotation,
    starts_with("sum_tot_n_snps.Female"), starts_with("sum_bi_n_snps.Female"),
    starts_with("sum_total_cells.Female"), starts_with("sum_cells_with_bi.Female"),
    starts_with("count_bi_individuals.Female"),starts_with("count_individuals.Female"),
    starts_with("sum_tot_n_snps.Male"), starts_with("sum_bi_n_snps.Male"),
    starts_with("sum_total_cells.Male"), starts_with("sum_cells_with_bi.Male"),
    starts_with("count_bi_individuals.Male"),starts_with("count_individuals.Male")
  )

summary_wide <- subset(summary_wide, !(sum_bi_n_snps.Female == 0 & sum_bi_n_snps.Male == 0))
summary_wide$check = ifelse(summary_wide$sum_cells_with_bi.Female > summary_wide$sum_cells_with_bi.Male, 'ok', 'check')
summary_wide$check[is.na(summary_wide$check)] = 'ok'
summary_wide = subset(summary_wide, !(check == 'check' & Domain.category != 'PAR'))
View(summary_wide)
summary_wide$key = paste(summary_wide$Gene, summary_wide$Annotation)

summary_df$key = paste(summary_df$Gene, summary_df$Annotation)
summary_df = subset(summary_df, key %in% summary_wide$key)
summary_df$sum_cells_with_bi[which(summary_df$sum_cells_with_bi == 0)] = NA
summary_df$Annotation = factor(summary_df$Annotation, levels = c(exn.levels, inn.levels, other.levels))
summary_df = subset(summary_df, Annotation %!in% c('Tcell','SM','Peri','L5ET','Fibro'))
genes_with_male_only <- summary_df %>%
  group_by(Gene) %>%
  filter(all(sex == "Male")) %>%
  pull(Gene)
summary_df <- summary_df %>%
  filter(!(Gene %in% genes_with_male_only & sex == "Male"))
summary_df$cellp = ifelse(summary_df$sex == 'Female',log(summary_df$sum_cells_with_bi+1),log(summary_df$sum_cells_with_bi+1)*-1)

summary_wide <- summary_wide %>%
  filter(!(Gene %in% genes_with_male_only))
summary_wide = subset(summary_wide, Annotation %!in% c('Tcell','SM','Peri','L5ET','Fibro'))

length(unique(summary_wide$Gene))
summary_wide %>% group_by(Domain.category) %>% summarise(length(unique(Gene)))
write.csv(summary_wide, file = 'summary_wide.csv')

ggplot(summary_df, aes(x = Gene, y = Annotation, fill = sum_cells_with_bi)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[2], 
                       high = sex.colors[1], 
                       na.value = 'white') +  
  facet_grid(sex~Domain.category, space = 'free', scales = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(summary_df, aes(x = Gene, y = Annotation, fill = cellp)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[2], 
                       high = sex.colors[1], 
                       limits = c(-5, 5),
                       na.value = 'white') +  
  facet_grid(sex~Domain.category, space = 'free', scales = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

xcie = data.frame(gene = unique(subset(summary_df, Domain.category!='PAR'))) 
sexres = readRDS('sexres.rds')
result <- xcie %>%
  rowwise() %>%
  mutate(
    ns = sum(sexres$gene == gene & is.na(sexres$bias), na.rm = T),
    female = sum(sexres$gene == gene & sexres$bias == 'female', na.rm = T),
    male = sum(sexres$gene == gene & sexres$bias == 'male', na.rm = T))
result <- result %>%
  rowwise() %>%
  mutate(max_category = names(.)[which.max(c(ns, female, male)) + 1])
table(result$max_category)
dim(subset(result, female > 0))

fnp = subset(summary_df, sex == 'Female' & Domain.category == 'nonPAR')
fnc = fnp %>% group_by(Gene) %>% summarise(n=n())
fnc = subset(fnc, n == 1)

summary_df = subset(summary_df, Gene %!in% fnc$Gene)

ggplot(summary_df, aes(x = Gene, y = Annotation, fill = cellp)) +
  geom_tile() +
  scale_fill_gradient2(low = sex.colors[2], 
                       high = sex.colors[1], 
                       na.value = 'white') +  
  facet_grid(sex~Domain.category, space = 'free', scales = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

head(summary_wide)

## comparison to previous

xcib = read.csv('xci-balaton.csv')
colnames(xcib)[1] = 'gene'
table(xcib$Domain.category)
table(xcib$Balaton.consensus.calls)
table(xcib$Cotton.AI)
table(xcib$Carrel.SNPs)
table(xcib$Cotton.DNAm)

xcib$gene[which(xcib$gene == 'KAL1')] = 'ANOS1'
xcib$gene[which(xcib$gene == 'BCLAF3')] = 'CXorf23'
xcib$gene[which(xcib$gene == 'GCNA')] = 'ACRC'
xcib$gene[which(xcib$gene == 'INTS6L')] = 'DDX26B'
xcib$gene[which(xcib$gene == 'NEXMIF')] = 'KIAA2022'

xcit = read.csv('xci-status.csv') #tukianen
colnames(xcit)[1] = 'gene'
xcib = merge(xcib, xcit[,c('gene','Combined.XCI.status')], by = 'gene', all=T)
table(xcib$Combined.XCI.status)

library(dplyr)
library(ggplot2)

xcie$gene[which(xcie$gene %!in% xcib$gene)]
xcib_filtered = merge(xcib, xcie, by = 'gene', all.y = T)
xcib_filtered = unique(xcib_filtered)
dim(xcib_filtered)

# Escape: meets all 3 criteria
escape_genes <- xcib_filtered %>%
  filter(
    Domain.category == "E" | Balaton.consensus.calls %in% c("E", "Mostly E") | Cotton.AI == "E" |
      Carrel.SNPs == 'E' | Cotton.DNAm %in% c('all F E','all F E or E+U','E or E+U and VE','E or E+U and VE and U') |
                                                Combined.XCI.status == 'escape'
  ) %>%
  distinct(gene) %>%
  pull(gene)
length(escape_genes)

# Variable escape: meets variable evidence criteria
variable_escape_genes <- xcib_filtered %>%
  filter(
    Balaton.consensus.calls %in% c("VE", "Mostly VE") | Cotton.AI == "VE" | Carrel.SNPs == 'VE' |
      Cotton.DNAm %in% c('all 4 states','S or S+U and VE','S or S+U and VE and U') | Combined.XCI.status == 'variable'
  ) %>%
  distinct(gene) %>%
  pull(gene)
variable_escape_genes = variable_escape_genes[which(variable_escape_genes %!in% escape_genes)]
length(variable_escape_genes)

all_genes <- unique(xcib_filtered$gene)

other_genes <- setdiff(all_genes, union(escape_genes, variable_escape_genes))

summary_df <- data.frame(
  Category = c("Escape", "Variable Escape", "Unreported"),
  Count = c(length(escape_genes), length(variable_escape_genes), length(other_genes))
)
summary_df

View(subset(xcib_filtered, gene %in% other_genes))
summary_df$Category = factor(summary_df$Category, levels = c('Unreported','Variable Escape','Escape'))
levels(summary_df$Category)

ggplot(summary_df, aes(y = "", x = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Escape" = "steelblue", 
                               "Variable Escape" = "orange", 
                               "Unreported" = "darkgray")) +
  theme_article() +
  theme(legend.position = 'top')

be = read.csv('HPA_brain_elevated.csv')
bs = read.csv('brain_specific_new.csv')

in_be <- other_genes %in% be$Gene
in_bs <- other_genes %in% bs$Gene
category <- ifelse(in_be, "in_be",'neither')
category <- ifelse(in_bs, "in_bs", category)
summary_df <- data.frame(Gene = other_genes, Category = category)
counts <- table(summary_df$Category)
counts
counts/sum(counts)*100
counts_df <- data.frame(
  Category = names(counts),
  Count = as.numeric(counts),
  stringsAsFactors = FALSE
)

f=subset(combined, Gene %in% other_genes)
View(f %>% group_by(Gene, sex) %>% summarise(sum = sum(cells_with_bi)))

ggplot(counts_df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = 'Set2')

