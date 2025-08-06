#!/bin/bash

# Load bedtools module
module load bedtools

# Path to genes.bed file (must be sorted and not contain chrM/chrY unless used)
GENES_BED_ORIG="genes.bed"
GENES_BED_FILTERED="genes.filtered.bed"

# Loop over all GWAS files
#for FILE in *.tsv.bgz; do
for FILE in *G43*.tsv.bgz; do
  echo "Processing $FILE..."

  BASE=$(basename "$FILE" .tsv.bgz)

  # Step 1: Filter significant variants and remove low confidence
  #zcat "$FILE" | awk -F'\t' 'NR==1 || ($12 < 5e-3 && $5 != "true")' > "${BASE}.sig.tsv"
  zcat "$FILE" | awk -F'\t' 'NR==1 || ($12 < 5e-05)' > "${BASE}.sig.tsv" # this output was used below****
  echo "  - Filtered significant variants"

  # Step 2: Convert to BED format
  awk -F'\t' 'NR>1 {
    split($1, a, ":");
    chrom = "chr"a[1];
    start = a[2] - 1;
    end = a[2];
    print chrom"\t"start"\t"end"\t"$1
  }' "${BASE}.sig.tsv" > "${BASE}.sig.bed"

  # Step 3: Get list of chromosomes
  awk -F'\t' 'NR>1 {
    split($1, a, ":");
    print "chr"a[1]
  }' "${BASE}.sig.tsv" | sort -u > "${BASE}.used_chromosomes.txt"

  # Step 4: Filter genes.bed for those chromosomes only
  grep -wFf "${BASE}.used_chromosomes.txt" "$GENES_BED_ORIG" > "$GENES_BED_FILTERED"

  # Step 5: Use bedtools closest to link variants to nearby genes
  bedtools closest -a "${BASE}.sig.bed" -b "$GENES_BED_FILTERED" -d > "${BASE}.closest.txt"

  # Step 6: Extract unique gene names from closest result
  cut -f8 "${BASE}.closest.txt" | sort | uniq > "${BASE}.closest.gene.txt"
  close_n=$(wc -l < "${BASE}.closest.gene.txt")
  echo "  - Found $close_n unique gene(s) from closest"

done
