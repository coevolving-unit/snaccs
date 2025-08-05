library(ggplot2)
library(stringr)
library(dplyr)
library(maditr)
library(egg)

#### load and clean data ####

gwas = read.csv('/GWAS/sex-biased-gwas.csv')
gwas = subset(gwas, Set != 'sex_biased')
unique(gwas$Disease)

gwas$code = str_sub(gwas$Disease, 1,3)
gwas$code2 = str_sub(gwas$Disease, -5,-1)
gwas = subset(gwas, code2 != '_both' & code2 != 'th_v2')
table(gwas$Disease)

#gwas = subset(gwas, p_value < 0.05)

gwas$code2 = as.factor(gwas$code2)
levels(gwas$code2) = c('male GWAS','female GWAS')

# only males for
# parkinsons_icd10_male
# autism_male
# remove
gwas = subset(gwas, code != 'aut')
gwas = subset(gwas, Disease != 'parkinsons_icd10_male')

gwas$Disease <- gsub("_female|_male", "", gwas$Disease)
gwas$Disease = as.factor(gwas$Disease)
levels(gwas$Disease) = c('AD (father)','AD (mother)','Mania','MS','OCD',
                         'PD (father)','PD (mother)',
                         'Schizophrenia')

#### end ####

### plot ####

# use ratio of ORs

po = reshape2::dcast(gwas, Cluster + Set + Disease ~ code2, value.var = c('OR')) 
pp = reshape2::dcast(gwas, Cluster + Set + Disease ~ code2, value.var = c('p_value')) 
po = cbind(po, pp$`male GWAS`, pp$`female GWAS`)
colnames(po) = c('cluster','set','disease','mgwasOR','fgwasOR','mgwasp','fgwasp')

po$cluster <- factor(po$cluster, levels = unique(po$cluster))
po$diff = log2(po$mgwasOR / po$fgwasOR)
po$code = ifelse(po$mgwasp < 0.05 | po$fgwasp < 0.05, 'keep', 'remove')
write.csv(po, file = 'gwas-results-ratio.csv')
po$colorcode = ifelse(po$set == 'female_biased' & po$diff < 0, 'female',
                      ifelse(po$set == 'male_biased' & po$diff > 0, 'male','ns'))
po = subset(po, code == 'keep')

table(po$diff > 0, po$mgwasp < 0.05)
table(po$diff < 0, po$fgwasp < 0.05)

ggplot(po, aes(y = disease, x = diff, color = colorcode))+
  facet_grid(cluster ~ set, drop = T) +
  scale_color_manual(values = c(sex.colors,"light grey")) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = 0, xend = diff, y = disease, yend = disease), size = 0.3) +
  geom_point(size = 2) +
  theme_article()

ggplot(po, aes(y = disease, x = diff, color = colorcode))+
  facet_grid(set ~ cluster, drop = T) +
  scale_color_manual(values = c(sex.colors,"light grey")) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = 0, xend = diff, y = disease, yend = disease), size = 0.3) +
  geom_point(size = 2) +
  theme_article() +
  theme(legend.position = 'bottom')

#### end ####
