#!/bin/bash
#SBATCH --job-name=process_bams
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00

# Process BAM files - extract X chromosome reads and rename barcodes

# Define directories
DIR="/data/cellranger_output/"
OUTDIR="/data/Xbams/"

# Load required modules
module load python/3.10 samtools

# Get the sample name from the batch array task ID
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: No sample name found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Define input/output BAM paths
INPUT_BAM="${DIR}${SAMPLE_NAME}/outs/possorted_genome_bam.bam"
OUTPUT_BAM="${OUTDIR}${SAMPLE_NAME}.barcode.renamed.X.bam"

# Ensure input BAM exists before proceeding
if [[ ! -f "$INPUT_BAM" ]]; then
    echo "ERROR: Input BAM file $INPUT_BAM not found!"
    exit 1
fi

echo "Processing $SAMPLE_NAME..."
echo "Input: $INPUT_BAM"
echo "Output: $OUTPUT_BAM"

# Run Python script to process BAM
python process_bam.py "$INPUT_BAM" "$OUTPUT_BAM"

# Index the output BAM
samtools index "$OUTPUT_BAM"

echo "Processing of ${SAMPLE_NAME} completed."
