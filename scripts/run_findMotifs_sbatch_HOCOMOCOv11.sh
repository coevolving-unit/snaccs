#!/bin/bash
#SBATCH --job-name=findMotifs_HOCOMOCOv11
#SBATCH --time=1:00:00
#SBATCH --mem=4GB

# Add HOMER to PATH
PATH=$PATH:/data/homer/bin/

# Run motif analysis using HOCOMOCOv11
findMotifs.pl $INPUT_FILE human $OUTPUT_PREFIX -start -1000 -end 300 -bg $BG_FILE -mset HOCOMOCOv11 -nomotif -nogo
