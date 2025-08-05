#!/bin/bash

# Add HOMER to PATH
PATH=$PATH:/data/homer/bin/
which findMotifs.pl

# List of clusters
clus=(1 2 3 4 5 6 7 8 9 10 11 12 13)

# Submit jobs via SLURM
for clus in "${clus[@]}"; do
    sbatch --export=CLUS=$clus,INPUT_FILE="cluster${clus}-homer-mn.txt",OUTPUT_PREFIX="cluster${clus}-homer-mn",BG_FILE="cluster${clus}-homer-bgn.txt" run_findMotifs_sbatch.sh
    sbatch --export=CLUS=$clus,INPUT_FILE="cluster${clus}-homer-fn.txt",OUTPUT_PREFIX="cluster${clus}-homer-fn",BG_FILE="cluster${clus}-homer-bgn.txt" run_findMotifs_sbatch.sh
    sbatch --export=CLUS=$clus,INPUT_FILE="cluster${clus}-homer-sn.txt",OUTPUT_PREFIX="cluster${clus}-homer-sn",BG_FILE="cluster${clus}-homer-bgn.txt" run_findMotifs_sbatch.sh
    sbatch --export=CLUS=$clus,INPUT_FILE="cluster${clus}-homer-mn.txt",OUTPUT_PREFIX="cluster${clus}-homer-mn-HOCOMOCO",BG_FILE="cluster${clus}-homer-bgn.txt" run_findMotifs_sbatch.sh
    sbatch --export=CLUS=$clus,INPUT_FILE="cluster${clus}-homer-fn.txt",OUTPUT_PREFIX="cluster${clus}-homer-fn-HOCOMOCO",BG_FILE="cluster${clus}-homer-bgn.txt" run_findMotifs_sbatch.sh
    sbatch --export=CLUS=$clus,INPUT_FILE="cluster${clus}-homer-sn.txt",OUTPUT_PREFIX="cluster${clus}-homer-sn-HOCOMOCO",BG_FILE="cluster${clus}-homer-bgn.txt" run_findMotifs_sbatch.sh
done
