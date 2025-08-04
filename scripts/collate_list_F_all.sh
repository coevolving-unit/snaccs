#!/bin/bash

# scripts/collate_list_F_all.sh  
# Collate results for female samples

#SBATCH --job-name=collate_females
#SBATCH --output=logs/collate_females_%A_%a.out
#SBATCH --error=logs/collate_females_%A_%a.err
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16

set -e

module load salmon/1.10.0

# Set up environment
export PATH="$(echo "$PATH" | sed -E 's!:*/usr/local/sbin:*!:!g; s/^://; s/:$//')"

# Activate conda environment if available
if command -v conda >/dev/null 2>&1; then
    source myconda 2>/dev/null || true
    mamba activate alevin-fry 2>/dev/null || conda activate alevin-fry 2>/dev/null || true
fi

# Get sample
sample_list=${1:-sample_list_F.txt}
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)

echo "Collating results for sample: $sample"

alevin-fry collate \
    -t 16 \
    -i mapped_reads_F/"${sample}"_quant \
    -r mapped_reads_F/"${sample}"_map

echo "Collation completed for sample: $sample"
