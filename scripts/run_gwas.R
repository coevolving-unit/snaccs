library(stringr)

sig = read.csv('/data/significant_genes_per_cluster.csv')
bg = read.csv('/data/background_genes_per_cluster.csv')
tx2gene = readRDS('/data/gene2chrom_mash.rds')

table(sig$bias, sig$cluster)

sig = subset(sig, gene %in% subset(tx2gene, chromosome_name %in% c(1:22))$external_gene_name)
bg = subset(bg, gene %in% subset(tx2gene, chromosome_name %in% c(1:22))$external_gene_name)

table(sig$bias, sig$cluster)
table(bg$cluster)

disease_files <- list.files(pattern = "*.closest.gene.txt")
disease_genes <- lapply(disease_files, function(f) unique(readLines(f)))
names(disease_genes) <- gsub("\\.closest\\.gene\\.txt", "", disease_files)

clus = unique(sig$cluster)

results <- data.frame(
  Cluster = character(),
  Set = character(),
  Disease = character(),
  OR = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:length(clus)){
  
  cl <- clus[i]
  cat("Processing cluster:", cl, "\n")
  
  mn = subset(sig, bias == 'male' & cluster == clus[i])$gene
  fn = subset(sig, bias == 'female' & cluster == clus[i])$gene
  bn = subset(sig, bias == 'both' & cluster == clus[i])$gene
  sn = unique(c(mn, fn, bn))
  bgn = subset(bg, cluster == clus[i])$gene
  
  gene_sets <- list(male_biased = mn,
                    female_biased = fn,
                    sex_biased = sn)
  
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    
    for (disease in names(disease_genes)) {
      disease_set <- disease_genes[[disease]]
      disease_set = disease_set[disease_set %in% bgn]
      
      # Perform Fisher's exact test
      overlap <- length(intersect(genes, disease_set))
      a <- overlap
      b <- length(genes) - a
      c <- length(disease_set) - a
      d <- length(setdiff(bgn, union(genes, disease_set)))  # background excluding both sets

      contingency <- matrix(c(a, b, c, d), nrow = 2)
      fisher <- fisher.test(contingency, alternative = 'greater')
      
  results <- rbind(results, data.frame(
        Cluster = cl,
        Set = set_name,
        Disease = disease,
        OR = fisher$estimate,
        p_value = fisher$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(results, 'sex-biased-gwas.csv')
                        
