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

/data/snRNAseq_fastq/
|-- {batch_id}/
|   |-- FASTQ/
|       |-- map_file.txt              
|       |-- {file_id}_R1_001.fastq.gz 
|       |-- {file_id}_R2_001.fastq.gz 
|       |-- ...
|-- {another_batch}/
|   |-- FASTQ/
|       |-- ...
  
# Pipeline
  
### Initial QC and mapping in CellRanger

* **Key libraries:** stringr

```
# run CellRanger
scripts/cellranger.sh
```

### Create Sex-Specific References

```
# Create sex-specific splici transcriptomes and references
bash scripts/sex-specific-splici-index.sh

# Build splici transcriptomes
Rscript scripts/make_splici_txome.R

# Build salmon indices
bash scripts/build_salmon_indices.sh
```

### Alevin-fry Quantification

```
# Create sample lists
# sample_list_M.txt - one male sample ID per line
# sample_list_F.txt - one female sample ID per line

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

```

### Create Count Matrices

```
# Create SingleCellExperiment objects
alevin_matrix_M.R
alevin_matrix_F.R

Rscript scripts/generate_matrix_swarms.R
swarm -f alevin_matrix_M_R.swarm --module R -g 10
swarm -f alevin_matrix_F_R.swarm --module R -g 10
```

### Quality Control and Filtering

```
# Generate QC reports
alevin_qc_m.R
alevin_qc_f.R

Rscript scripts/generate_matrix_swarms.R
swarm -f alevin_qc_m.swarm -g 50 --module R
swarm -f alevin_qc_f.swarm -g 50 --module R
```

### Merge Metadata and Further Processing

```
# Merge sample metadata
Rscript scripts/merge_metadata.R

# Remove ambient RNA using CellBender
bash scripts/cell-bender.sh

# Remove doublets and low-quality cells
Rscript scripts/mito-doublets.R
```

### Integration and Clustering

```
# Iterative integration, clustering, and filtering
Rscript scripts/integration-clean.R
```

### Cell Type Annotation

```
# Prepare reference datasets
bash scripts/allen_M1_ref.sh
bash scripts/allen_MCA_ref.sh  
bash scripts/allen_MTG_ref.sh
bash scripts/BRAIN_ref.sh

# Prepare data and run annotation
bash scripts/prep-for-annotation.sh
bash scripts/annotation-swarm.sh

# Visualize annotation results
Rscript scripts/plot_azimuth.R
```

### Pseudobulk Analysis

```
# Create pseudobulk samples for downstream analysis
Rscript scripts/pseudobulk.R
```

### Proportions

* **Key libraries:** stringr

```
# Test for sex differences in cell type proportions
Rscript scripts/plot_calc_proportions.R
```

### Differential Expression Analysis

```
# Run differential expression using limma-voom
bash scripts/differential_expression.sh

# Run multivariate adaptive shrinkage (MASHR)
bash scripts/mashr.sh
```

### Transcriptome-wide Impact Analysis

```
# Estimate sex-specific transcriptome-wide impact (sexTWI)
bash scripts/TRADE.sh

# Leave-one-out validation
bash scripts/DE-loo-TRADE-jk.sh

# Visualize results
Rscript scripts/TRADE.R
```

### Variance Partitioning

```
# Partition variance across biological and technical factors
bash scripts/variance-partitioning.sh
Rscript scripts/plot-variance-partitioning.R
```

### Functional Enrichment Analysis

```
# Cluster autosomal sex effects
Rscript scripts/plot-sex-effects.R

# Transcription factor enrichment using HOMER
bash scripts/homer-clusters.sh

# Gene set enrichment analysis (GSEA)
bash scripts/gsea-analysis.sh

# GWAS enrichment analysis
bash scripts/gwas-enrichment.sh
```

### Allele-Specific Expression

```
# Estimate allele-specific expression patterns
bash scripts/ase-analysis.sh
Rscript scripts/plot-ase-results.R
```

### Output Structure

results/
├── cellranger/              # Cell Ranger outputs
├── alevin/                  # Alevin-fry quantification
│   ├── mapped_reads_M/      # Male samples
│   └── mapped_reads_F/      # Female samples
├── qc/                      # Quality control reports
├── integration/             # Integrated data objects
├── annotation/              # Cell type annotations
├── pseudobulk/             # Pseudobulk expression matrices
├── differential_expression/ # DE analysis results
├── functional_enrichment/   # Pathway and TF enrichments
├── variance_partitioning/   # Variance decomposition results
└── figures/                # Publication-ready plots

