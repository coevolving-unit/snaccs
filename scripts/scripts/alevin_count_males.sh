#!/bin/bash

# scripts/alevin_count_males.sh
# Process male samples with Alevin-fry

#SBATCH --job-name=alevin_males
#SBATCH --output=logs/alevin_males_%A_%a.out
#SBATCH --error=logs/alevin_males_%A_%a.err
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16

set -e

module load salmon

# Set paths
index=af_tutorial_splici/grch38_splici_idx
sample_list=${1:-sample_list_M.txt}

# Get current sample ID
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)
echo "Processing sample: ${sample}"

# Create sample directory
mkdir -p mapped_reads_M/"${sample}"

# Function to find FASTQ files for sample
# Customize this function based on your file organization
find_fastq_files() {
    local sample_id=$1
    local fastq_base_dir=${2:-/path/to/fastq/files}
    
    # Method 1: Direct pattern matching (customize paths)
    find $fastq_base_dir -name "*${sample_id}*R1*.fastq.gz" | sort > mapped_reads_M/"${sample}"/"${sample}"_R1.txt
    find $fastq_base_dir -name "*${sample_id}*R2*.fastq.gz" | sort > mapped_reads_M/"${sample}"/"${sample}"_R2.txt
    
    # Method 2: Using metadata file (uncomment and customize if you have one)
    # if [[ -f "sample_metadata.txt" ]]; then
    #     batch=$(awk -v sample="$sample_id" '$1==sample {print $2}' sample_metadata.txt)
    #     find ${fastq_base_dir}/${batch} -name "*${sample_id}*R1*.fastq.gz" | sort > mapped_reads_M/"${sample}"/"${sample}"_R1.txt
    #     find ${fastq_base_dir}/${batch} -name "*${sample_id}*R2*.fastq.gz" | sort > mapped_reads_M/"${sample}"/"${sample}"_R2.txt
    # fi
}

# Find FASTQ files (customize the base directory path)
FASTQ_BASE_DIR="/path/to/your/fastq/files"
find_fastq_files "$sample" "$FASTQ_BASE_DIR"

# Convert file lists to space-delimited format
r1_files=$(cat mapped_reads_M/"${sample}"/"${sample}"_R1.txt | tr '\n' ' ')
r2_files=$(cat mapped_reads_M/"${sample}"/"${sample}"_R2.txt | tr '\n' ' ')

# Check that files were found
if [[ -z "$r1_files" || -z "$r2_files" ]]; then
    echo "Error: No FASTQ files found for sample $sample"
    exit 1
fi

echo "Found R1 files: $r1_files"
echo "Found R2 files: $r2_files"

# Run Alevin
echo "Running Alevin mapping for sample: $sample"
salmon alevin \
    -i $index \
    -p 16 \
    -l IU \
    --chromiumV3 \
    --sketch \
    -1 $r1_files \
    -2 $r2_files \
    -o mapped_reads_M/"${sample}"_map \
    --tgMap transcriptome_splici_fl86/transcriptome_splici_fl86_t2g.tsv

echo "Alevin mapping completed for sample: $sample"
