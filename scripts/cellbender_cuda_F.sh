#!/bin/bash

# scripts/cellbender_cuda_F.sh
# Run CellBender on female samples

#SBATCH --job-name=cellbender_F
#SBATCH --output=logs/cellbender_F_%A_%a.out
#SBATCH --error=logs/cellbender_F_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k80:1
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16

set -e

# Get sample ID
sample_list=${1:-sample_list_F.txt}
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list)

echo "Running CellBender on female sample: $sample"

# Set up working directory
WORKDIR="mapped_reads_F/${sample}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || { echo "Failed to enter directory $WORKDIR"; exit 1; }

# Load CellBender module
module load cellbender

# Check if input file exists
input_file="../${sample}.h5ad"
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file not found: $input_file"
    exit 1
fi

echo "Input file: $input_file"
echo "Output file: ${sample}_bender.h5"

# Run CellBender
cellbender remove-background \
    --cuda \
    --input "$input_file" \
    --output "${sample}_bender.h5" \
    --expected-cells 5000 \
    --total-droplets 20000 \
    --checkpoint-mins 60 \
    --learning-rate 0.000025 \
    --epochs 200 \
    --fpr 0

echo "CellBender completed for sample: $sample"
