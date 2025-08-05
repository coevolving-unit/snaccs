#!/bin/bash

# scripts/cellranger.sh
# Run Cell Ranger count on samples

#SBATCH --job-name=cellranger
#SBATCH --output=logs/cellranger_%A_%a.out
#SBATCH --error=logs/cellranger_%A_%a.err
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8

set -e

module load cellranger/5.0.1

# Determine sex and sample from array task ID
# First get total counts for each sex
MALE_COUNT=0
FEMALE_COUNT=0

if [[ -f "sample_list_M.txt" ]]; then
    MALE_COUNT=$(wc -l < sample_list_M.txt)
fi

if [[ -f "sample_list_F.txt" ]]; then
    FEMALE_COUNT=$(wc -l < sample_list_F.txt)
fi

# Determine which sample list to use based on array task ID
if (( SLURM_ARRAY_TASK_ID <= MALE_COUNT )); then
    # Process male sample
    sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list_M.txt)
    TRANSCRIPTOME="/path/to/your/GRCh38_male_reference"  # Male reference with Y
    SEX="M"
    echo "Processing male sample: $sample"
    echo "Using male reference (with Y chromosome)"
elif (( SLURM_ARRAY_TASK_ID <= MALE_COUNT + FEMALE_COUNT )); then
    # Process female sample
    FEMALE_TASK_ID=$((SLURM_ARRAY_TASK_ID - MALE_COUNT))
    sample=$(sed -n ${FEMALE_TASK_ID}p sample_list_F.txt)
    TRANSCRIPTOME="/path/to/your/GRCh38_noY_reference"   # Female reference without Y
    SEX="F"
    echo "Processing female sample: $sample"
    echo "Using female reference (Y-masked)"
else
    echo "Error: Array task ID $SLURM_ARRAY_TASK_ID exceeds total sample count"
    exit 1
fi

# Function to find FASTQ directories and sample names using metadata approach
find_fastq_info_metadata() {
    local sample_id=$1
    local metadata_file=${2:-"metadata/region_samples_list.txt"}
    local fastq_base_dir=${3:-"/data/snRNAseq_fastq"}
    
    echo "Finding FASTQ info for sample: $sample_id using metadata approach"
    
    # Create temporary directory for this sample
    temp_dir="cellranger_temp_${sample_id}"
    mkdir -p "$temp_dir"
    
    # Get batch information for current sample
    if [[ ! -f "$metadata_file" ]]; then
        echo "Error: Metadata file not found: $metadata_file"
        echo "Please ensure region_samples_list.txt is available"
        exit 1
    fi
    
    # Extract batch info (assuming column 13 is sample_id and column 10 is batch)
    awk -F'\t' -v sample="$sample_id" 'NR == 1 || $13 == sample' "$metadata_file" > "${temp_dir}/${sample_id}_batch.txt"
    
    if [[ ! -s "${temp_dir}/${sample_id}_batch.txt" ]]; then
        echo "Error: No batch information found for sample $sample_id"
        exit 1
    fi
    
    # Get batch ID
    awk -F'\t' '{print $10}' "${temp_dir}/${sample_id}_batch.txt" | sort | uniq > "${temp_dir}/${sample_id}_batch2.txt"
    sed -i '1d' "${temp_dir}/${sample_id}_batch2.txt"  # Remove header
    
    batch=$(cat "${temp_dir}/${sample_id}_batch2.txt")
    
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
    
    awk -v sample="$sample_id" '$1 == sample' "$map_file" > "${temp_dir}/${sample_id}.txt"
    
    if [[ ! -s "${temp_dir}/${sample_id}.txt" ]]; then
        echo "Error: No entries found in map file for sample $sample_id"
        exit 1
    fi
    
    # Get unique file identifiers (these become sample names for Cell Ranger)
    awk '{print $2}' "${temp_dir}/${sample_id}.txt" | sort | uniq > "${temp_dir}/${sample_id}_uniq.txt"
    
    # Set FASTQ directory and sample names
    FASTQ_DIRS="${fastq_base_dir}/${batch}/FASTQ"
    SAMPLE_NAMES=$(cat "${temp_dir}/${sample_id}_uniq.txt" | tr '\n' ',' | sed 's/,$//')
    
    echo "FASTQ directory: $FASTQ_DIRS"
    echo "Sample names: $SAMPLE_NAMES"
    
    # Verify that FASTQ files exist for at least one sample name
    first_sample_name=$(echo "$SAMPLE_NAMES" | cut -d',' -f1)
    if ! ls "${FASTQ_DIRS}"/*"${first_sample_name}"* >/dev/null 2>&1; then
        echo "Error: No FASTQ files found for sample name: $first_sample_name in $FASTQ_DIRS"
        exit 1
    fi
    
    # Clean up temporary directory
    rm -rf "$temp_dir"
}

# Create output directory
mkdir -p cellranger_output

# Run Cell Ranger count
echo "Running Cell Ranger count for sample: $sample"

cellranger count \
    --id="${sample}" \
    --transcriptome="$TRANSCRIPTOME" \
    --fastqs="$FASTQ_DIRS" \
    --sample="$SAMPLE_NAMES" \
    --chemistry=SC3Pv3 \
    --include-introns \
    --localcores=8 \
    --localmem=60

# Move output to organized directory
if [[ -d "${sample}" ]]; then
    mv "${sample}" cellranger_output/
    echo "Cell Ranger output moved to cellranger_output/${sample}"
fi

echo "Cell Ranger processing completed for sample: $sample"
