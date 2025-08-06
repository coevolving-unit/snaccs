#!/bin/bash

# scripts/permit_list_M_all.sh
# Generate permit lists for male samples

#SBATCH --job-name=permit_males
#SBATCH --output=logs/permit_males_%A_%a.out
#SBATCH --error=logs/permit_males_%A_%a.err
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4

set -e

module load salmon/1.10.0

# Set up environment
export PATH="$(echo "$PATH" | sed -E 's!:*/usr/local/sbin:*!:!g; s/^://; s/:$//')"

# Activate conda environment if available
if command -v conda >/dev/null 2>&1; then
    source myconda 2>/dev/null || true
    mamba activate alevin-fry 2>/dev/null || conda activate alevin-fry 2>/dev/null || true
fi

# Get sample list
sample_list=${1:-sample_list_M.txt}
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)

echo "Generating permit list for sample: $sample"

# Create output directory
mkdir -p results/alevin/permit_lists

# Download barcode whitelist if not present in data directory
mkdir -p data/references/barcodes
if [[ ! -f "data/references/barcodes/3M-february-2018.txt" ]]; then
    echo "Downloading 10x barcode whitelist..."
    wget https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz -O data/references/barcodes/3M-february-2018.txt.gz
    gunzip data/references/barcodes/3M-february-2018.txt.gz
fi

# Generate permit list
echo "Running alevin-fry generate-permit-list..."
alevin-fry generate-permit-list \
    -d fw \
    -i results/alevin/mapped_reads_M/"${sample}"_map \
    -o results/alevin/mapped_reads_M/"${sample}"_quant \
    --unfiltered-pl data/references/barcodes/3M-february-2018.txt

echo "Permit list generation completed for sample: $sample"
echo "Output saved to: results/alevin/mapped_reads_M/${sample}_quant"
