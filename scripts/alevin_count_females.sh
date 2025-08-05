#!/bin/bash

# scripts/alevin_count_females.sh
# Process female samples with Alevin-fry

#SBATCH --job-name=alevin_females
#SBATCH --output=logs/alevin_females_%A_%a.out
#SBATCH --error=logs/alevin_females_%A_%a.err
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16

set -e

module load salmon

# Set paths
index=af_tutorial_splici_noY/grch38_splici_idx
sample_list=${1:-sample_list_F.txt}

# Get current sample ID
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)
echo "Processing sample: ${sample}"

# Create sample directory
mkdir -p mapped_reads_F/"${sample}"

# Function to find FASTQ files using metadata approach
find_fastq_files_metadata() {
    local sample_id=$1
    local metadata_file=${2:-"metadata/region_samples_list.txt"}
    local fastq_base_dir=${3:-"/data/snRNAseq_fastq"}
    
    echo "Finding FASTQ files for sample: $sample_id using metadata approach"
    
    # Get batch information for current sample
    if [[ ! -f "$metadata_file" ]]; then
        echo "Error: Metadata file not found: $metadata_file"
        echo "Please ensure region_samples_list.txt is available"
        exit 1
    fi
    
    # Extract batch info (assuming column 13 is sample_id and column 10 is batch)
    awk -F'\t' -v sample="$sample_id" 'NR == 1 || $13 == sample' "$metadata_file" > mapped_reads_F/"${sample}"_batch.txt
    
    if [[ ! -s mapped_reads_F/"${sample}"_batch.txt ]]; then
        echo "Error: No batch information found for sample $sample_id"
        exit 1
    fi
    
    # Get batch ID
    awk -F'\t' '{print $10}' mapped_reads_F/"${sample}"_batch.txt | sort | uniq > mapped_reads_F/"${sample}"_batch2.txt
    sed -i '1d' mapped_reads_F/"${sample}"_batch2.txt  # Remove header
    
    batch=$(cat mapped_reads_F/"${sample}"_batch2.txt)
    
    if [[ -z "$batch" ]]; then
        echo "Error: Could not determine batch for sample $sample_id"
        exit 1
    fi
    
    echo "Found batch: $batch"
    
    # Get files for sample using map file
    map_file="${fastq_base_dir}/${batch}/FASTQ/map_file.txt"
    if [[ ! -f "$map_file" ]]; then
        echo "Error: Map file not found: $map_file"
        exit 1
    fi
    
    awk -v sample="$sample_id" '$1 == sample' "$map_file" > mapped_reads_F/"${sample}".txt
    
    if [[ ! -s mapped_reads_F/"${sample}".txt ]]; then
        echo "Error: No entries found in map file for sample $sample_id"
        exit 1
    fi
    
    # Get unique file identifiers
    awk '{print $2}' mapped_reads_F/"${sample}".txt | sort | uniq > mapped_reads_F/"${sample}"_uniq.txt
    
    # Find R1 and R2 files
    fastq_dir="${fastq_base_dir}/${batch}/FASTQ"
    
    for f in $(cat mapped_reads_F/"${sample}"_uniq.txt); do 
        ls "$fastq_dir" | grep -E "${f}\." | grep -E "R1" || true
    done > mapped_reads_F/"${sample}"_R1.txt
    
    for f in $(cat mapped_reads_F/"${sample}"_uniq.txt); do 
        ls "$fastq_dir" | grep -E "${f}\." | grep -E "R2" || true
    done > mapped_reads_F/"${sample}"_R2.txt
    
    # Check if files were found
    if [[ ! -s mapped_reads_F/"${sample}"_R1.txt || ! -s mapped_reads_F/"${sample}"_R2.txt ]]; then
        echo "Error: No R1/R2 FASTQ files found for sample $sample_id"
        exit 1
    fi
    
    # Add full path prefix
    prefix="${fastq_dir}/"
    awk -v prefix="$prefix" '{print prefix $0}' mapped_reads_F/"${sample}"_R1.txt > mapped_reads_F/"${sample}"_R1_prefix.txt
    awk -v prefix="$prefix" '{print prefix $0}' mapped_reads_F/"${sample}"_R2.txt > mapped_reads_F/"${sample}"_R2_prefix.txt
    
    # Convert to space-delimited format
    cat mapped_reads_F/"${sample}"_R1_prefix.txt | xargs -rd'\n' > mapped_reads_F/"${sample}"_R1_space.txt
    cat mapped_reads_F/"${sample}"_R2_prefix.txt | xargs -rd'\n' > mapped_reads_F/"${sample}"_R2_space.txt
    
    # Verify files exist
    r1_files=$(cat mapped_reads_F/"${sample}"_R1_space.txt)
    r2_files=$(cat mapped_reads_F/"${sample}"_R2_space.txt)
    
    echo "Found R1 files: $r1_files"
    echo "Found R2 files: $r2_files"
    
    # Check that at least one file exists for each
    first_r1=$(echo $r1_files | awk '{print $1}')
    first_r2=$(echo $r2_files | awk '{print $1}')
    
    if [[ ! -f "$first_r1" || ! -f "$first_r2" ]]; then
        echo "Error: FASTQ files not found at expected paths"
        echo "First R1 file: $first_r1"
        echo "First R2 file: $first_r2"
        exit 1
    fi
}

# Run Alevin
echo "Running Alevin mapping for sample: $sample"
salmon alevin \
    -i $index \
    -p 16 \
    -l IU \
    --chromiumV3 \
    --sketch \
    -1 $(cat mapped_reads_F/"${sample}"_R1_space.txt) \
    -2 $(cat mapped_reads_F/"${sample}"_R2_space.txt) \
    -o mapped_reads_F/"${sample}"_map \
    --tgMap transcriptome_splici_fl_noY86/transcriptome_splici_fl86_t2g.tsv

echo "Alevin mapping completed for sample: $sample"
