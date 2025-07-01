library(stringr)
library(dplyr)

# load metadata
# brain-level metrics (PMI, ph)
# sample-level metrics (RIN) 

meta_data = read.delim('/data/region_samples_list.txt')
meta_data$RIN = as.numeric(meta_data$RIN)
meta_data$orig.ident = str_sub(meta_data$sample_id, -5, -1)
meta_data$new.id = str_sub(meta_data$new_id, -5, -1)

# sequencing_metrics.csv is from cellranger output
# sample-level metrics

seq_info = read.csv('/data/sequencing_metrics.csv')
seq_info = subset(seq_info, type == 'all')
meta_data = merge(meta_data, seq_info, by = 'new.id')
colnames(meta_data)[2] = 'orig.ident'
meta_data$orig.ident = str_sub(meta_data$orig.ident, -5, -1)
meta_data$Number.of.Reads = gsub(",","",meta_data$Number.of.Reads)
meta_data$Number.of.Reads = as.numeric(meta_data$Number.of.Reads)
meta_data$Reads.Mapped.to.Genome = gsub("%","",meta_data$Reads.Mapped.to.Genome)
meta_data$Reads.Mapped.to.Genome = as.numeric(meta_data$Reads.Mapped.to.Genome)
meta_data$Estimated.Number.of.Cells = gsub(",","",meta_data$Estimated.Number.of.Cells)
meta_data$Estimated.Number.of.Cells = as.numeric(meta_data$Estimated.Number.of.Cells)

# alevin_qc_tables.csv is from AlevinQC
# sample-level metrics (mapping rate, # reads)

seq_info = read.csv('/data/alevin_qc_tables.csv')
s = data.frame(orig.ident = seq_info$sample,
				nreads = seq_info$Total.number.of.processed.reads,
				maprate = seq_info$Number.of.mapped.reads / seq_info$Total.number.of.processed.reads)
s$orig.ident = str_sub(s$orig.ident, -5, -1)

# get cell counts from cleaned/filtered data
# sample-level metrics (total_cells)

other_annotate_class_counts = read.csv('outputs/other_annotate_class_counts_corrected.csv')
inn_annotate_class_counts = read.csv('outputs/inn_annotate_class_counts_corrected.csv')
exn_upper_annotate_class_counts = read.csv('outputs/exn_upper_annotate_class_counts_corrected.csv')
exn_lower_annotate_class_counts = read.csv('outputs/exn_lower_annotate_class_counts_corrected.csv')

all_counts = cbind(other_annotate_class_counts,
				inn_annotate_class_counts[,-1],
				exn_upper_annotate_class_counts[,-1],
				exn_lower_annotate_class_counts[,-1])
all_counts$total_cells = rowSums(all_counts[,c(2:length(colnames(all_counts)))])
colnames(all_counts)[1] = 'new.id'

# merge metadata
# average/sum sample-level stats for merged samples 
# RIN is same for replicates
# total cells are already estimated using new_id
# (average = mapping rate; sum = # reads)

meta_data = meta_data[,-1]
meta_all = merge(meta_data, s, by = 'orig.ident', all = T)
meta_all = merge(meta_all, all_counts, by = 'new.id', all = T)
dim(meta_all)

m = data.frame(meta_all %>% group_by(new.id) %>% summarise(n=n()))
d = subset(m, n == 2)
r = subset(meta_all, new.id %in% d$new.id)
n = data.frame(r %>% group_by(new.id) %>% summarise(nreads = sum(nreads), maprate = mean(maprate)))

for(i in 1:length(n$new.id)){
	meta_all[which(meta_all$new.id == n$new.id[i]),c('nreads','maprate')] = n[i,c(2:3)]
	}
meta_all[which(meta_all$orig.ident == 79557),'batch'] = meta_all[which(meta_all$orig.ident == 81708),'batch']

saveRDS(meta_all, file = 'outputs/meta_merged.rds')

