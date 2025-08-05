#!/usr/bin/env Rscript

library(stringr)

# create 1 file for each cell type x region (columns = samples)

pseudo = read.csv('other_class_pseudo_corrected.csv') # done
pseudo = read.csv('inn_class_pseudo_corrected.csv') # done
pseudo = read.csv('exn_upper_class_pseudo_corrected.csv') # done
pseudo = read.csv('exn_lower_class_pseudo_corrected.csv') # done

meta = read.csv('/data/region_samples_list.csv')
meta$new_id = str_sub(meta$new_id, -5, -1)
regions = unique(meta$region)

rownames(pseudo) = pseudo$X
pseudo = pseudo[,-1]

# class
n = unique(gsub("^.+?\\.(.+?)\\..*$", "\\1", colnames(pseudo)))
n
for(i in 1:length(n)){
  print(n[i])
  for(j in 1:length(regions)){
    dnow = pseudo[,which(gsub("^.+?\\.(.+?)\\..*$", "\\1", colnames(pseudo)) == n[i])]
    s = subset(meta, region == regions[j])$new_id
    dnow = dnow[,which(sub('.*\\.', '', colnames(dnow)) %in% s)]
    colnames(dnow) = sub('.*\\.', '', colnames(dnow))
    saveRDS(dnow, file = paste(regions[j],n[i],'pseudo_corrected.rds',sep="-"))
  }
}

# create 1 file containing all cell types x regions x samples (columns = cell types x regions x samples)
# this is used for overall variance partitioning (exp ~ class + region + individual + sex + etc)

c1 = read.csv('other_class_pseudo_corrected.csv')
c2 = read.csv('inn_class_pseudo_corrected.csv')
c3 = read.csv('exn_upper_class_pseudo_corrected.csv')
c4 = read.csv('exn_lower_class_pseudo_corrected.csv')

rownames(c1) = c1$X
c1 = c1[,-1]
rownames(c2) = c2$X
c2 = c2[,-1]
rownames(c3) = c3$X
c3 = c3[,-1]
rownames(c4) = c4$X
c4 = c4[,-1]

# combine exn upper and lower
a = cbind(c3, c4)

saveRDS(a, 'exn_class_pseudo_corrected.rds') # done

# combine all
a = cbind(c1, c2, c3, c4)

saveRDS(a, 'all_class_pseudo_corrected.rds') # done

# create 1 file containing all sample-level pseudo (use class!)

a = readRDS('all_class_pseudo_corrected.rds')

s = unique(sub('.*\\.', '', colnames(a)))
s
o = matrix(0, nrow = 36591, ncol = length(s))
rownames(o) = rownames(a)
colnames(o) = s

length(s)
for(i in 1:length(s)) {
  print(i)
  dnow = a[,which(sub('.*\\.', '', colnames(a)) == s[i])]
  su = rowSums(dnow)
  o[,i] = su
}

saveRDS(o, 'all_sample_level_pseudo_corrected.rds') # done
