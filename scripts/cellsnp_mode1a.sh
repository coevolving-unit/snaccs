#!/bin/bash
#SBATCH --job-name=cellsnp_1a
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --time=24:00:00

# Run cellsnp-lite in mode 1a (with cell barcodes, using known SNPs from mode 2b)

module load cellsnp-lite

# Get sample names
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)
SAMPLE_NAME_SHORT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams_short.txt)

if [ -z "$SAMPLE_NAME" ] || [ -z "$SAMPLE_NAME_SHORT" ]; then
    echo "ERROR: Sample name not found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Define paths
BARCODE_RNA="barcodes/${SAMPLE_NAME}_barcodes.txt"
BAM_RNA="Xbams/${SAMPLE_NAME}.barcode.renamed.X.bam"
VCF_REF="Xcalls/${SAMPLE_NAME}.RNA.chrX.2b/cellSNP.cells.vcf"
OUTPUT_DIR="Xcalls/sample_${SAMPLE_NAME_SHORT}.RNA.chrX.1a"
REFSEQ="chromosome.X.fa"
UMITAG_RNA="UB"

echo "Running cellsnp-lite mode 1a for $SAMPLE_NAME..."
echo "Input BAM: $BAM_RNA"
echo "Barcodes: $BARCODE_RNA"
echo "Reference VCF: $VCF_REF"
echo "Output directory: $OUTPUT_DIR"

# Check if required input files exist
for file in "$BAM_RNA" "$BARCODE_RNA" "$VCF_REF" "$REFSEQ"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Check if output already exists
if [ -d "$OUTPUT_DIR" ] && [ -f "$OUTPUT_DIR/cellSNP.base.vcf" ]; then
    echo "Output already exists for $SAMPLE_NAME, skipping..."
    exit 0
fi

# Run cellsnp-lite mode 1a
cellsnp-lite \
    -s "$BAM_RNA" \
    -b "$BARCODE_RNA" \
    --cellTAG CB \
    --UMItag "$UMITAG_RNA" \
    -p 10 \
    -R "$VCF_REF" \
    --refseq "$REFSEQ" \
    --chrom=chrX \
    --minCOUNT 5 \
    --minMAF 0.05 \
    --minMAPQ 5 \
    -O "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    echo "cellsnp-lite mode 1a completed successfully for $SAMPLE_NAME"
    echo "Results saved to: $OUTPUT_DIR"
    
    # Report some basic statistics
    if [ -f "$OUTPUT_DIR/cellSNP.base.vcf" ]; then
        VCF_LINES=$(wc -l < "$OUTPUT_DIR/cellSNP.base.vcf")
        echo "Generated VCF with $VCF_LINES lines"
    fi
else
    echo "ERROR: cellsnp-lite mode 1a failed for $SAMPLE_NAME"
    exit 1
fi
