library(stringr)

samples_list = read.csv('/data/region_samples_list.csv')
regions = levels(as.factor(samples_list$region))

class_short = "other"	
subtype_now = c('Oligo','OPC','Astro','Micro','BEC')

for (i in 1:length(regions)) {
	for(k in 1:length(subtype_now)) {
		line = paste("Rscript voom_limma.R", regions[i], subtype_now[k], sep = " ")
		write(line, file = "limma_other_class.swarm", append = TRUE, sep = "\n")
	}}
	
class_short = "inn"	
subtype_now = c('ADARB2','Chandelier','VIP','SST','LAMP5','PVALB','PAX6')

for (i in 1:length(regions)) {
	for(k in 1:length(subtype_now)) {
		line = paste("Rscript voom_limma.R", regions[i], subtype_now[k], sep = " ")
		write(line, file = "limma_inn_class.swarm", append = TRUE, sep = "\n")
	}}

class_short = "exc"
subtype_now = c('L23IT','L4IT','L5IT','L6IT','L56NP','L6b','L6CT')

for (i in 1:length(regions)) {
	for(k in 1:length(subtype_now)) {
		line = paste("Rscript voom_limma.R", regions[i], subtype_now[k], sep = " ")
		write(line, file = "limma_exn_class.swarm", append = TRUE, sep = "\n")
	}}
