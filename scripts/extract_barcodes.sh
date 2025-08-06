#!/bin/bash
#SBATCH --job-name=extract_barcodes
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=6:00:00

# Extract cell barcodes from CellBender H5 files
# Handles swapped sample names by using bams_short.txt for correct mapping

module load hdf5

# Get sample names
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)
SAMPLE_NAME_SHORT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams_short.txt)

if [ -z "$SAMPLE_NAME" ] || [ -z "$SAMPLE_NAME_SHORT" ]; then
    echo "ERROR: Sample name not found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Extracting barcodes for $SAMPLE_NAME (using $SAMPLE_NAME_SHORT)"

# Define paths
H5_FILE="/data/*/$SAMPLE_NAME_SHORT/${SAMPLE_NAME_SHORT}_bender_filtered.h5"
BARCODE_TMP="barcodes/${SAMPLE_NAME}_barcodes_tmp.txt"
BARCODE_FINAL="barcodes/${SAMPLE_NAME}_barcodes.txt"

# Create barcodes directory
mkdir -p barcodes

# Check if final barcode file already exists
if [ -f "$BARCODE_FINAL" ]; then
    echo "Barcode file already exists: $BARCODE_FINAL"
    exit 0
fi

# Extract barcodes from H5 file
echo "Extracting barcodes from H5 file..."
h5dump -d /matrix/barcodes -o "$BARCODE_TMP" $H5_FILE

if [ ! -f "$BARCODE_TMP" ]; then
    echo "ERROR: Failed to extract barcodes from H5 file"
    exit 1
fi

# Process the extracted barcodes
echo "Processing barcodes..."
sed 's/.*: //; s/, /\n/g; s/"//g; s/,$//' "$BARCODE_TMP" > "$BARCODE_FINAL"

# Remove header line and add -RNA suffix
sed -i '1d' "$BARCODE_FINAL"
sed -i 's/$/-RNA/' "$BARCODE_FINAL"

# Clean up temporary file
rm -f "$BARCODE_TMP"

# Report results
BARCODE_COUNT=$(wc -l < "$BARCODE_FINAL")
echo "Extracted $BARCODE_COUNT barcodes to: $BARCODE_FINAL"

echo "Barcode extraction completed for $SAMPLE_NAME"
