#############
# download annotations
#############

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

gunzip gencode.v44.annotation.gtf.gz

awk '$3 == "gene"' gencode.v44.annotation.gtf | \
awk 'BEGIN{OFS="\t"} {
  match($0, /gene_name "([^"]+)"/, arr); 
  print $1, $4-1, $5, arr[1]
}' > genes.bed

grep -v -w "chrM" genes.bed > genes.nochrM.bed

grep -v -w "chrY" genes.nochrM.bed > genes.nochrM.noY.bed

awk '$3 == "gene"' gencode.v44.annotation.gtf | \
awk 'BEGIN{OFS="\t"} {
  match($0, /gene_id "([^"]+)"/, id);
  match($0, /gene_name "([^"]+)"/, name);
  print id[1], name[1]
}' > ensembl_to_symbol.tsv

sort ensembl_to_symbol.tsv > sorted_mapping.tsv

#############
# download UKBB GWAS results
#############

# Illnesses of father: Alzheimer's disease/dementia
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20107_10.gwas.imputed_v3.both_sexes.tsv.bgz -O 20107_10.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20107_10.gwas.imputed_v3.female.tsv.bgz -O 20107_10.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20107_10.gwas.imputed_v3.male.tsv.bgz -O 20107_10.gwas.imputed_v3.male.tsv.bgz

# Illnesses of mother: Alzheimer's disease/dementia
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20110_10.gwas.imputed_v3.both_sexes.tsv.bgz -O 20110_10.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20110_10.gwas.imputed_v3.female.tsv.bgz -O 20110_10.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20110_10.gwas.imputed_v3.male.tsv.bgz -O 20110_10.gwas.imputed_v3.male.tsv.bgz

# Illnesses of father: Parkinson's disease
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20107_11.gwas.imputed_v3.both_sexes.tsv.bgz -O 20107_11.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20107_11.gwas.imputed_v3.female.tsv.bgz -O 20107_11.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20107_11.gwas.imputed_v3.male.tsv.bgz -O 20107_11.gwas.imputed_v3.male.tsv.bgz

# Illnesses of mother: Parkinson's disease
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20110_11.gwas.imputed_v3.both_sexes.tsv.bgz -O 20110_11.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20110_11.gwas.imputed_v3.female.tsv.bgz -O 20110_11.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20110_11.gwas.imputed_v3.male.tsv.bgz -O 20110_11.gwas.imputed_v3.male.tsv.bgz

# Diagnoses - main ICD10: G20 Parkinson's disease
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G20.gwas.imputed_v3.both_sexes.tsv.bgz -O G20.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G20.gwas.imputed_v3.male.tsv.bgz -O G20.gwas.imputed_v3.male.tsv.bgz

# Mental health problems ever diagnosed by a professional: Autism, Asperger's or autistic spectrum disorder
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_14.gwas.imputed_v3.both_sexes.tsv.bgz -O 20544_14.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_14.gwas.imputed_v3.male.tsv.bgz -O 20544_14.gwas.imputed_v3.male.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_14.gwas.imputed_v3.both_sexes.v2.tsv.bgz -O 20544_14.gwas.imputed_v3.both_sexes.v2.tsv.bgz

# Mental health problems ever diagnosed by a professional: Obsessive compulsive disorder (OCD)
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_7.gwas.imputed_v3.both_sexes.tsv.bgz -O 20544_7.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_7.gwas.imputed_v3.female.tsv.bgz -O 20544_7.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_7.gwas.imputed_v3.male.tsv.bgz -O 20544_7.gwas.imputed_v3.male.tsv.bgz

# Schizophrenia, schizotypal and delusional disorders
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/F5_SCHIZO.gwas.imputed_v3.both_sexes.tsv.bgz -O F5_SCHIZO.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/F5_SCHIZO.gwas.imputed_v3.female.tsv.bgz -O F5_SCHIZO.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/F5_SCHIZO.gwas.imputed_v3.male.tsv.bgz -O F5_SCHIZO.gwas.imputed_v3.male.tsv.bgz

# Mental health problems ever diagnosed by a professional: Mania, hypomania, bipolar or manic-depression
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_10.gwas.imputed_v3.both_sexes.v2.tsv.bgz -O 20544_10.gwas.imputed_v3.both_sexes.v2.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_10.gwas.imputed_v3.female.v2.tsv.bgz -O 20544_10.gwas.imputed_v3.female.v2.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20544_10.gwas.imputed_v3.male.v2.tsv.bgz -O 20544_10.gwas.imputed_v3.male.v2.tsv.bgz

# Diagnoses - main ICD10: G35 Multiple sclerosis
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G35.gwas.imputed_v3.both_sexes.tsv.bgz -O G35.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G35.gwas.imputed_v3.female.tsv.bgz -O G35.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G35.gwas.imputed_v3.male.tsv.bgz -O G35.gwas.imputed_v3.male.tsv.bgz

# Diagnoses - main ICD10: G43 Migraine
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G43.gwas.imputed_v3.both_sexes.tsv.bgz -O G43.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G43.gwas.imputed_v3.female.tsv.bgz -O G43.gwas.imputed_v3.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/G43.gwas.imputed_v3.male.tsv.bgz -O G43.gwas.imputed_v3.male.tsv.bgz
