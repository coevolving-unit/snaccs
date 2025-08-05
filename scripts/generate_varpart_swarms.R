library(stringr)

samples_list = read.csv('/dataregion_samples_list.csv')
regions = levels(as.factor(samples_list$region))

class_short = "other"	
subtype_now = c('Oligo','OPC','Astro','Fibro','Micro','BEC','SM','Peri','Tcell')

for (i in 1:length(regions)) {
	for(k in 1:length(subtype_now)) {
		line = paste("Rscript varpart.R", regions[i], subtype_now[k], sep = " ")
		write(line, file = "varpart_other_class.swarm", append = TRUE, sep = "\n")
	}}
	
class_short = "inn"	
subtype_now = c('ADARB2','Chandelier','VIP','SST','LAMP5','PVALB','PAX6')

for (i in 1:length(regions)) {
	for(k in 1:length(subtype_now)) {
		line = paste("Rscript varpart.R", regions[i], subtype_now[k], sep = " ")
		write(line, file = "varpart_inn_class.swarm", append = TRUE, sep = "\n")
	}}

class_short = "exc"
subtype_now = c('L23IT','L4IT','L5IT','L6IT','L5ET','L56NP','L6b','L6CT')

for (i in 1:length(regions)) {
	for(k in 1:length(subtype_now)) {
		line = paste("Rscript varpart.R", regions[i], subtype_now[k], sep = " ")
		write(line, file = "varpart_exn_class.swarm", append = TRUE, sep = "\n")
	}}
