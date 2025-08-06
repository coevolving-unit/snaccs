#!/bin/bash

# scripts/sex-specific-splici-index.sh
# Creates sex-specific splici transcriptome indices

#SBATCH --job-name=splici_index
#SBATCH --output=logs/splici_index_%j.out
#SBATCH --error=logs/splici_index_%j.err
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8

set -e

echo "Creating sex-specific splici transcriptome indices..."

# Create standardized directory structure
mkdir -p data/references data/references/source results/references

# Download 10x reference to data directory
cd data/references
if [[ ! -f "refdata-gex-GRCh38-2020-A.tar.gz" ]]; then
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
fi
if [[ ! -d "refdata-gex-GRCh38-2020-A" ]]; then
    tar xzf refdata-gex-GRCh38-2020-A.tar.gz
fi

# Genome metadata
genome="GRCh38_noY"
version="2020-A"

# Set up directories
build="GRCh38-2020-A_build_noY"
source="source"
mkdir -p "$build" "$source"

cd "$source"

# Download reference files
fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
if [[ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]]; then
    echo "Downloading primary assembly FASTA..."
    wget $fasta_url
    gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa
fi

gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
if [[ ! -f "gencode.v32.primary_assembly.annotation.gtf" ]]; then
    echo "Downloading GENCODE GTF..."
    wget $gtf_url
    gunzip -c gencode.v32.primary_assembly.annotation.gtf.gz > gencode.v32.primary_assembly.annotation.gtf
fi

cd ..

fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"

# Modify FASTA headers to match GENCODE format
echo "Modifying FASTA headers..."
fasta_modified="$build/$(basename "$fasta_in").modified"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffixes from GTF IDs
echo "Processing GTF annotations..."
gtf_modified="$build/$(basename "$gtf_in").modified"
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define biotype patterns
BIOTYPE_PATTERN="(protein_coding|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""

# Create gene allowlist
echo "Creating gene allowlist..."
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

echo "Gene allowlist contains $(wc -l < ${build}/gene_allowlist) genes"

# Filter GTF
echo "Filtering GTF..."
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$gtf_modified" >> "$gtf_filtered"

# Create Y chromosome mask
echo "Creating Y chromosome mask..."
module load bioawk bedtools

bioawk -c fastx '$name=="chrY" {print $name"\t1\t"length($seq)}' "$fasta_modified" > rm_Y.bed

# Mask Y chromosome for female reference
echo "Masking Y chromosome for female reference..."
fasta_modified_noY="$build/$(basename "$fasta_in").modified.maskedY"
bedtools maskfasta -fi "$fasta_modified" -bed rm_Y.bed -fo "$fasta_modified_noY"

# Build Cell Ranger references
echo "Building Cell Ranger references..."
module load cellranger/6.0.0

# Male reference (with Y) - output to data/references
cellranger mkref \
    --ref-version="$version" \
    --genome="GRCh38_male" \
    --fasta="$fasta_modified" \
    --genes="$gtf_filtered"

# Female reference (Y masked) - output to data/references  
cellranger mkref \
    --ref-version="$version" \
    --genome="$genome" \
    --fasta="$fasta_modified_noY" \
    --genes="$gtf_filtered"

# Unzip GTF files
echo "Unzipping GTF files..."
gunzip GRCh38_male/genes/genes.gtf.gz
gunzip GRCh38_noY/genes/genes.gtf.gz

echo "Sex-specific references created successfully in data/references/"
echo "Male reference: data/references/GRCh38_male"  
echo "Female reference: data/references/GRCh38_noY"
