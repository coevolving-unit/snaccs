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

# Configuration - use environment variables or defaults
METADATA_FILE=${REGION_SAMPLES_LIST:-"data/metadata/region_samples_list.txt"}
FASTQ_BASE_DIR=${FASTQ_BASE_PATH:-"data/fastq"}
REFERENCE_DIR=${REFERENCE_PATH:-"data/references"}

# Set paths using standardized structure
index="${REFERENCE_DIR}/af_tutorial_splici/grch38_splici_idx"
sample_list=${1:-sample_list_M.txt}

# Get current sample ID
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)
echo "Processing sample: ${sample}"

# Create sample directory in standardized location
mkdir -p results/alevin/mapped_reads_M/"${sample}"

# Function to find FASTQ files using metadata approach
find_fastq_files_metadata() {
    local sample_id=$1
    local metadata_file=$2
    local fastq_base_dir=$3
    
    echo "Finding FASTQ files for sample: $sample_id using metadata approach"
    
    # Get batch information for current sample
    if [[ ! -f "$metadata_file" ]]; then
        echo "Error: Metadata file not found: $metadata_file"
        echo "Please ensure region_samples_list.txt is available in data/metadata/"
        exit 1
    fi
    
    # Extract batch info (assuming column 13 is sample_id and column 10 is batch)
    awk -F'\t' -v sample="$sample_id" 'NR == 1 || $13 == sample' "$metadata_file" > results/alevin/mapped_reads_M/"${sample}"_batch.txt
    
    if [[ ! -s results/alevin/mapped_reads_M/"${sample}"_batch.txt ]]; then
        echo "Error: No batch information found for sample $sample_id"
        exit 1
    fi
    
    # Get batch ID
    awk -F'\t' '{print $10}' results/alevin/mapped_reads_M/"${sample}"_batch.txt | sort | uniq > results/alevin/mapped_reads_M/"${sample}"_batch2.txt
    sed -i '1d' results/alevin/mapped_reads_M/"${sample}"_batch2.txt  # Remove header
    
    batch=$(cat results/alevin/mapped_reads_M/"${sample}"_batch2.txt)
    
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
    
    awk -v sample="$sample_id" '$1 == sample' "$map_file" > results/alevin/mapped_reads_M/"${sample}".txt
    
    if [[ ! -s results/alevin/mapped_reads_M/"${sample}".txt ]]; then
        echo "Error: No entries found in map file for sample $sample_id"
        exit 1
    fi
    
    # Get unique file identifiers
    awk '{print $2}' results/alevin/mapped_reads_M/"${sample}".txt | sort | uniq > results/alevin/mapped_reads_M/"${sample}"_uniq.txt
    
    # Find R1 and R2 files
    fastq_dir="${fastq_base_dir}/${batch}/FASTQ"
    
    for f in $(cat results/alevin/mapped_reads_M/"${sample}"_uniq.txt); do 
        ls "$fastq_dir" | grep -E "${f}\." | grep -E "R1" || true
    done > results/alevin/mapped_reads_M/"${sample}"_R1.txt
    
    for f in $(cat results/alevin/mapped_reads_M/"${sample}"_uniq.txt); do 
        ls "$fastq_dir" | grep -E "${f}\." | grep -E "R2" || true
    done > results/alevin/mapped_reads_M/"${sample}"_R2.txt
    
    # Check if files were found
    if [[ ! -s results/alevin/mapped_reads_M/"${sample}"_R1.txt || ! -s results/alevin/mapped_reads_M/"${sample}"_R2.txt ]]; then
        echo "Error: No R1/R2 FASTQ files found for sample $sample_id"
        exit 1
    fi
    
    # Add full path prefix
    prefix="${fastq_dir}/"
    awk -v prefix="$prefix" '{print prefix $0}' results/alevin/mapped_reads_M/"${sample}"_R1.txt > results/alevin/mapped_reads_M/"${sample}"_R1_prefix.txt
    awk -v prefix="$prefix" '{print prefix $0}' results/alevin/mapped_reads_M/"${sample}"_R2.txt > results/alevin/mapped_reads_M/"${sample}"_R2_prefix.txt
    
    # Convert to space-delimited format
    cat results/alevin/mapped_reads_M/"${sample}"_R1_prefix.txt | xargs -rd'\n' > results/alevin/mapped_reads_M/"${sample}"_R1_space.txt
    cat results/alevin/mapped_reads_M/"${sample}"_R2_prefix.txt | xargs -rd'\n' > results/alevin/mapped_reads_M/"${sample}"_R2_space.txt
    
    # Verify files exist
    r1_files=$(cat results/alevin/mapped_reads_M/"${sample}"_R1_space.txt)
    r2_files=$(cat results/alevin/mapped_reads_M/"${sample}"_R2_space.txt)
    
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

# Find FASTQ files for this sample
find_fastq_files_metadata "$sample" "$METADATA_FILE" "$FASTQ_BASE_DIR"

# Run Alevin
echo "Running Alevin mapping for sample: $sample"
salmon alevin \
    -i "$index" \
    -p 16 \
    -l IU \
    --chromiumV3 \
    --sketch \
    -1 $(cat results/alevin/mapped_reads_M/"${sample}"_R1_space.txt) \
    -2 $(cat results/alevin/mapped_reads_M/"${sample}"_R2_space.txt) \
    -o results/alevin/mapped_reads_M/"${sample}"_map \
    --tgMap "${REFERENCE_DIR}/transcriptome_splici_fl86/transcriptome_splici_fl86_t2g.tsv"

echo "Alevin mapping completed for sample: $sample"
