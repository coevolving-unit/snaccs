#!/bin/bash

# scripts/permit_list_F_all.sh
# Generate permit lists for female samples

#SBATCH --job-name=permit_females
#SBATCH --output=logs/permit_females_%A_%a.out
#SBATCH --error=logs/permit_females_%A_%a.err
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
sample_list=${1:-sample_list_F.txt}
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)

echo "Generating permit list for sample: $sample"

# Download barcode whitelist if not present
if [[ ! -f "3M-february-2018.txt" ]]; then
    wget https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz
    gunzip 3M-february-2018.txt.gz
fi

# Generate permit list
alevin-fry generate-permit-list \
    -d fw \
    -i mapped_reads_F/"${sample}"_map \
    -o mapped_reads_F/"${sample}"_quant \
    --unfiltered-pl 3M-february-2018.txt

echo "Permit list generation completed for sample: $sample"
