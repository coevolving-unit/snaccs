tx2gene = readRDS('gene2chrom_mash.rds')

sig = read.csv('significant_genes_per_cluster.csv')
bg = read.csv('background_genes_per_cluster.csv')

# autosomes only
sig = subset(sig, gene %in% subset(tx2gene, chromosome_name %in% c(1:22))$external_gene_name)
bg = subset(bg, gene %in% subset(tx2gene, chromosome_name %in% c(1:22))$external_gene_name)

clus = unique(sig$cluster)

for(i in 1:length(clus)){
  
  print(clus[i])
  mn = subset(sig, bias == 'male' & cluster == clus[i])$gene
  fn = subset(sig, bias == 'female' & cluster == clus[i])$gene
  bn = subset(sig, bias == 'both' & cluster == clus[i])$gene
  sn = unique(c(mn, fn, bn))
  bgn = subset(bg, cluster == clus[i])$gene
  
  write.table(fn, file = paste('cluster',clus[i],'-homer-fn.txt',sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(mn, file = paste('cluster',clus[i],'-homer-mn.txt',sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(sn, file = paste('cluster',clus[i],'-homer-sn.txt',sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(bgn, file = paste('cluster',clus[i],'-homer-bgn.txt',sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE)
 
  }
