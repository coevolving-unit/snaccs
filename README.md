# snaccs

#### Repository for our Single Nuclei Analysis of Cerebral Cortex Sex differences (SNACCS)

This repository contains scripts used in the analysis of sex effects on snRNAseq data from human brains,

Note that we ran most steps on the NIH ([Biowulf](https://hpc.nih.gov/)) high-performance computing cluster. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Biowulf uses the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.3). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

* Raw FASTQ files: 10x Chromium single-nuclei RNA-seq data
* Sample metadata: Tab-delimited file with sample information including sex
* Sample lists:
* sample_list_M.txt: List of XY sample IDs
* sample_list_F.txt: List of XX sample IDs
* Reference data: Human genome (GRCh38) and annotation files

### Directory Structure

```
project_root/
├── data/
│   ├── fastq/                        # Raw FASTQ files
│   ├── references/                   # Reference genomes and indices
│   └── metadata/                     # Sample metadata files
├── results/
│   ├── cellranger/                   # Cell Ranger outputs
│   ├── alevin/                       # Alevin-fry quantification
│   ├── cellbender/                   # Ambient RNA removal
│   ├── qc/                           # Quality control reports
│   ├── integration/                  # Integrated datasets
│   ├── annotation/                   # Cell type annotations
│   ├── pseudobulk/                   # Pseudobulk expression
│   ├── differential_expression/       # DE analysis results
│   ├── variance_partitioning/        # Variance analysis
│   ├── functional_enrichment/        # GO/pathway analysis
│   ├── gwas_enrichment/              # GWAS enrichment
│   └── allelic_expression/           # X chromosome analysis
├── scripts/                          # All analysis scripts
└── logs/                             # Job log files
```

### FASTQ Directory Structure

```
data/fastq/
├── {batch_id}/
│   └── FASTQ/
│       ├── map_file.txt              # Maps sample IDs to file identifiers
│       ├── {file_id}_R1_001.fastq.gz # Read 1 files
│       ├── {file_id}_R2_001.fastq.gz # Read 2 files
│       └── ...
├── {another_batch}/
│   └── FASTQ/
│       └── ...
```
  
# Pipeline

### Create Sex-Specific References

```
# Create sex-specific splici transcriptomes and references
bash scripts/sex-specific-splici-index.sh

# Build splici transcriptomes
Rscript scripts/make_splici_txome.R

# Build salmon indices
bash scripts/build_salmon_indices.sh
```

### Create sample lists
```
# sample_list_M.txt - one male sample ID per line
# sample_list_F.txt - one female sample ID per line
```

### Cell Ranger QC and quantification

```
# Run Cell Ranger count
sbatch --array=1-$(wc -l < sample_list_M.txt) --mem=64G --time=8:00:00 --cpus-per-task=8 scripts/cellranger_males.sh
sbatch --array=1-$(wc -l < sample_list_F.txt) --mem=64G --time=8:00:00 --cpus-per-task=8 scripts/cellranger_females.sh
```

### Alevin-fry Quantification

```
# Map reads using Alevin-fry
sbatch --array=1-$(wc -l sample_list_M.txt | cut -d ' ' -f 1) --mem=300g --time=12:00:00 scripts/alevin_count_males.sh
sbatch --array=1-$(wc -l sample_list_F.txt | cut -d ' ' -f 1) --mem=300g --time=12:00:00 scripts/alevin_count_females.sh

# Generate permit lists
sbatch --array=1-$(wc -l sample_list_M.txt | cut -d ' ' -f 1) scripts/permit_list_M_all.sh
sbatch --array=1-$(wc -l sample_list_F.txt | cut -d ' ' -f 1) scripts/permit_list_F_all.sh

# Collate results
sbatch --array=1-$(wc -l sample_list_M.txt | cut -d ' ' -f 1) scripts/collate_list_M_all.sh
sbatch --array=1-$(wc -l sample_list_F.txt | cut -d ' ' -f 1) scripts/collate_list_F_all.sh

# Quantify gene expression
sbatch --array=1-$(wc -l sample_list_M.txt | cut -d ' ' -f 1) scripts/quantify_M.sh
sbatch --array=1-$(wc -l sample_list_F.txt | cut -d ' ' -f 1) scripts/quantify_F.sh

# Create SingleCellExperiment objects
# Individual scripts: alevin_matrix_M.R, alevin_matrix_F.R
Rscript scripts/generate_matrix_swarms.R
swarm -f alevin_matrix_M_R.swarm --module R -g 10
swarm -f alevin_matrix_F_R.swarm --module R -g 10

# Generate QC reports
# Individual scripts: alevin_qc_m.R, alevin_qc_f.R
Rscript scripts/generate_matrix_swarms.R
swarm -f alevin_qc_m.swarm -g 50 --module R
swarm -f alevin_qc_f.swarm -g 50 --module R
```

### Ambient RNA 

```
# Convert SCE objects to h5ad format for CellBender
# Individual scripts: sce-to-h5ad-M.R, sce-to-h5ad-F.R
Rscript scripts/generate_conversion_swarms.R
swarm -f sce-to-h5ad-M.swarm --module R/4.3 -g 10
swarm -f sce-to-h5ad-F.swarm --module R/4.3 -g 10

# Run CellBender for ambient RNA removal (requires GPU)
sbatch --array=1-$(wc -l sample_list_M.txt | cut -d ' ' -f 1) scripts/cellbender_cuda_M.sh
sbatch --array=1-$(wc -l sample_list_F.txt | cut -d ' ' -f 1) scripts/cellbender_cuda_F.sh
```

### Merge Metadata 

```
# Merge sample metadata
Rscript scripts/merge_metadata.R
```

### QC
```
# Convert CellBender results to Seurat objects
# Individual scripts: bender-seu-M.R, bender-seu-F.R
Rscript scripts/generate_bender_swarms.R
swarm -f bender-seu-M.swarm -g 200 --module R/4.3
swarm -f bender-seu-F.swarm -g 200 --module R/4.3

# Remove mitochondrial cells and doublets, generate QC plots
# Individual scripts: qc-M.R, qc-F.R
Rscript scripts/generate_qc_swarms.R
swarm -f qc-M.swarm -g 15 --module R/4.3
swarm -f qc-F.swarm -g 15 --module R/4.3
```

### Integration and Clustering

```
# Step 1: Convert Seurat objects to AnnData format
# Individual scripts: seurat-to-anndata_F.R; seurat-to-anndata_M.R
Rscript scripts/generate_anndata_swarms.R
swarm -f seurat-to-anndata_M.swarm -g 15 --module R/4.3
swarm -f seurat-to-anndata_F.swarm -g 15 --module R/4.3

# Step 2: Initial data integration and clustering
python scripts/03_data_integration.py

# Step 3: Cell type annotation based on marker genes  
python scripts/04_cell_type_annotation.py

# Step 4: Detailed subclustering of specific cell types
# Use template for custom cell types:
python scripts/05_subclustering_template.py <input_file> <cell_type> <output_prefix>

# Step 5: Excitatory neuron analysis (4 rounds of iterative clustering)
python scripts/06_excitatory_neurons.py

# Step 6: Inhibitory neuron analysis (3 rounds of iterative clustering)  
python scripts/07_inhibitory_neurons.py

# Step 7: Non-neuronal cell analysis (4 rounds of iterative clustering)
python scripts/08_nonneuronal_cells.py

# Step 8: Final integration and visualization
python scripts/09_final_integration.py
```

### Cell Type Annotation

```
# Download referendes
bash scripts/references.sh

# Prepare references
Rscript scripts/prepare_references.R

# Prepare data 
python scripts/prep-for-annotation.py

# Generate swarms and run annotation
# Individual file: azimuth.R
Rscript scripts/generate_azimuth_swarm.R
swarm -f azimuth.swarm --module R/4.3 -g 200

# Visualize annotation results
Rscript scripts/plot_azimuth.R
```

### Pseudobulk Analysis

```
# Create pseudobulk samples for downstream analysis
python scripts/pseudobulk.py
Rscript scripts/pseudobulk.R
```

### Proportions
```
# Test for sex differences in cell type proportions
Rscript scripts/plot_calc_proportions.R
```

### Differential Expression Analysis

```
# Run differential expression using limma-voom
Rscript scripts/voom_limma.R
Rscript scripts/generate_voom_swarms.R
swarm -f limma_other_class.swarm -g 20 --module R/4.3 
swarm -f limma_inn_class.swarm -g 20 --module R/4.3  
swarm -f limma_exn_class.swarm -g 20 --module R/4.3  

# Process DE results
Rscript scripts/process_de_results.R

# Run multivariate adaptive shrinkage (MASHR)
Rscript scripts/mashr_sex_class_approach1.R
Rscript scripts/mashr_sex_class_approach2.R
Rscript scripts/mashr_sex_class_approach3.R
Rscript scripts/mashr_sex_class_approach4.R
Rscript scripts/mashr_sex_class_approach5.R
```

### Transcriptome-wide Impact Analysis

```
# Estimate sex-specific transcriptome-wide impact (sexTWI)
Rscript scripts/TRADE.R
Rscript scripts/TRADE-leave-one-out.R
Rscript scripts/TRADE-perm.R
Rscript scripts/generate_TRADE_swarms.R
swarm -f TRADE_other_class.swarm --module R/4.3 
swarm -f TRADE_inn_class.swarm --module R/4.3 
swarm -f TRADE_exn_class.swarm --module R/4.3
swarm -f TRADE-other-class-leave-one-out.swarm --module R/4.3 
swarm -f TRADE-inn-class-leave-one-out.swarm --module R/4.3 
swarm -f TRADE-exn-class-leave-one-out.swarm --module R/4.3 
swarm -f TRADE-other-class-perm.swarm --module R/4.3 
swarm -f TRADE-inn-class-perm.swarm --module R/4.3 
swarm -f TRADE-exn-class-perm.swarm --module R/4.3 

# Process and combine results
Rscript scripts/process_TRADE.R

# Visualize results
Rscript scripts/plot_TRADE.R
```

### Variance Partitioning

```
# Partition variance within cell types
Rscript scripts/varpart.R
Rscript scripts/generate_varpart_swarms.R
swarm -f varpart_exn_class.swarm --module R/4.3 
swarm -f varpart_inn_class.swarm --module R/4.3 
swarm -f varpart_other_class.swarm --module R/4.3 

# Process results
Rscript scripts/process_varpart.R

# Partition variance overall
Rscript scripts/varpart_all.R

# Plot results
Rscript scripts/plot-varpart.R
```

### Plot sex effects

```
# Plot sex effects + cluster autosomal sex effects
Rscript scripts/plot-sex-effects.R
```

### TF Enrichment Analysis (autosomal clusters)

```
# Get gene lists for HOMER
Rscript scripts/get_homer_gene_lists.R

# Run HOMER
bash scripts/run_findMotifs_vertebrates.sh
bash scripts/run_findMotifs_HOCOMOCOv11.sh
bash scripts/run_findMotifs_sbatch_vertebrates.sh
bash scripts/run_findMotifs_sbatch_HOCOMOCOv11.sh

# Plot HOMER results
Rscript scripts/plot_homer.R
```

### Functional Enrichment Analysis (autosomal clusters)
```
# Gene ontology (GO) enrichment analysis
Rscript scripts/run_plot_GO.R
```

### Sex-specific GWAS Enrichment Analysis (autosomal clusters)
```
# Download sex-specific GWAS 
bash scripts/download_gwas.sh

# Process sex-specific GWAS
bash scripts/process_gwas.sh
chmod +x process_gwas.sh
./process_gwas.sh

# Run GWAS enrichments
Rscript scripts/run_gwas.R

# Visualize sex-specific GWAS enrichment
Rscript scripts/plot_gwas.R
```

### Allele-Specific Expression

```
# Estimate allele-specific expression
bash scripts/setup_genotyping.sh
bash scripts/process_bams.sh
python scripts/process_bams.py
bash scripts/cellsnp_mode2b.sh
bash scripts/extract_barcodes.sh
bash scripts/cellsnp_mode1a.sh
Rscript scripts/process_cellsnp_results.R
bash scripts/run_annovar.sh
Rscript scripts/qc_analysis.R
Rscript scripts/combine_results.R
Rscript scripts/plot_allelic_escape.R
```
