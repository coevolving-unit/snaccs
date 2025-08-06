#!/bin/bash
#SBATCH --job-name=annovar
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=6:00:00

# Run ANNOVAR annotation on processed cellsnp results

module load annovar

# Get sample name
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: No sample name found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Define paths
INPUT_TSV="Xcalls/${SAMPLE_NAME}.RNA.chrX.1a.summary.tsv.gz"
ANNOVAR_INPUT="annovar/${SAMPLE_NAME}.RNA.chrX.snp.call.summary.annov.input"
ANNOVAR_OUTPUT="annovar/${SAMPLE_NAME}.RNA.chrX.snp.call"

echo "Running ANNOVAR for $SAMPLE_NAME..."

# Check if input file exists
if [[ ! -f "$INPUT_TSV" ]]; then
    echo "ERROR: Input TSV file not found: $INPUT_TSV"
    echo "Please run 06_process_cellsnp_results.R first"
    exit 1
fi

# Check if output already exists
if [[ -f "${ANNOVAR_OUTPUT}.hg38_multianno.txt" ]]; then
    echo "ANNOVAR output already exists for $SAMPLE_NAME, skipping..."
    exit 0
fi

# Create annovar directory
mkdir -p annovar

echo "Preparing ANNOVAR input from $INPUT_TSV..."

# Convert TSV to ANNOVAR input format
# Format: chr start end ref alt variant_id
zcat "$INPUT_TSV" | \
awk -F"\t" 'NR>1 {print $2"\t"$3"\t"$3"\t"$4"\t"$5"\t"$1}' | \
LANG=C sort | LANG=C uniq > "$ANNOVAR_INPUT"

INPUT_VARIANTS=$(wc -l < "$ANNOVAR_INPUT")
echo "Prepared $INPUT_VARIANTS unique variants for annotation"

if [ $INPUT_VARIANTS -eq 0 ]; then
    echo "ERROR: No variants found in input file"
    exit 1
fi

echo "Running ANNOVAR annotation..."

# Run ANNOVAR
table_annovar.pl \
    "$ANNOVAR_INPUT" \
    "$ANNOVAR_DATA/hg38" \
    -buildver hg38 \
    -out "$ANNOVAR_OUTPUT" \
    -remove \
    -protocol refGene,ensGene,avsnp150,ALL.sites.2015_08 \
    -operation g,g,f,f \
    -nastring NA \
    -polish

if [ $? -eq 0 ]; then
    echo "ANNOVAR completed successfully for $SAMPLE_NAME"
    
    # Report basic statistics
    if [[ -f "${ANNOVAR_OUTPUT}.hg38_multianno.txt" ]]; then
        ANNOTATED_VARIANTS=$(tail -n +2 "${ANNOVAR_OUTPUT}.hg38_multianno.txt" | wc -l)
        echo "Annotated $ANNOTATED_VARIANTS variants"
    fi
    
    # Clean up intermediate files
    rm -f "$ANNOVAR_INPUT"
    rm -f "${ANNOVAR_OUTPUT}.avinput"
    rm -f "${ANNOVAR_OUTPUT}.hg38_ALL.sites.2015_08_dropped"
    rm -f "${ANNOVAR_OUTPUT}.hg38_avsnp150_dropped"
else
    echo "ERROR: ANNOVAR failed for $SAMPLE_NAME"
    exit 1
fi

echo "ANNOVAR annotation completed for $SAMPLE_NAME"
