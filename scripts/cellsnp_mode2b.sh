#!/bin/bash
#SBATCH --job-name=cellsnp_2b
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=12:00:00

# Run cellsnp-lite in mode 2b (without cell barcodes, for genotyping)

module load cellsnp-lite

# Get sample name
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: No sample name found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Define paths
BAM_RNA="Xbams/${SAMPLE_NAME}.barcode.renamed.X.bam"
OUTPUT_DIR="Xcalls/${SAMPLE_NAME}.RNA.chrX.2b"
UMITAG_RNA="UB"

# Check if input BAM exists
if [[ ! -f "$BAM_RNA" ]]; then
    echo "ERROR: Input BAM file $BAM_RNA not found!"
    exit 1
fi

echo "Running cellsnp-lite mode 2b for $SAMPLE_NAME..."
echo "Input BAM: $BAM_RNA"
echo "Output directory: $OUTPUT_DIR"

# Run cellsnp-lite mode 2b
cellsnp-lite \
    -s "$BAM_RNA" \
    --minMAF 0.05 \
    --minCOUNT 5 \
    --cellTAG None \
    --UMItag "$UMITAG_RNA" \
    -p 10 \
    --genotype \
    --chrom=chrX \
    -O "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    echo "cellsnp-lite mode 2b completed successfully for $SAMPLE_NAME"
    echo "Results saved to: $OUTPUT_DIR"
else
    echo "ERROR: cellsnp-lite mode 2b failed for $SAMPLE_NAME"
    exit 1
fi
