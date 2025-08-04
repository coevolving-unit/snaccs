#!/bin/bash

# scripts/build_salmon_indices.sh
# Build salmon indices for both male and female splici transcriptomes

#SBATCH --job-name=salmon_index
#SBATCH --output=logs/salmon_index_%j.out
#SBATCH --error=logs/salmon_index_%j.err
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16

set -e

module load salmon

echo "Building salmon indices for splici transcriptomes"

# Create directories
mkdir -p af_tutorial_splici af_tutorial_splici_noY

# Build male index
echo "Building male (with Y chromosome) salmon index"
cd af_tutorial_splici
salmon index \
    -t ../transcriptome_splici_fl86/transcriptome_splici_fl86.fa \
    -i grch38_splici_idx \
    -p 16

cd ..

# Build female index  
echo "Building female (Y-masked) salmon index"
cd af_tutorial_splici_noY
salmon index \
    -t ../transcriptome_splici_fl_noY86/transcriptome_splici_fl86.fa \
    -i grch38_splici_idx \
    -p 16

cd ..

echo "Salmon index building completed successfully"
