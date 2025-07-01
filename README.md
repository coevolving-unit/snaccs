# snaccs

#### Repository for our Single Nuclei Analysis of Cerebral Cortex Sex differences (SNACCS)

This repository contains scripts used in the analysis of sex effects on snRNAseq data from human brains,

Note that we ran most steps on the NIH ([Biowulf](https://hpc.nih.gov/)) high-performance computing cluster. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Biowulf uses the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.3). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

* XXX
  
# Pipeline
  
### initial QC and mapping in CellRanger

* **Key libraries:** stringr

```
# run CellRanger
scripts/cellranger.sh
```

### Make splici transcriptomes

* **Key libraries:** stringr

```
# make sex specific splici transcriptomes
scripts/make_splici_txome.sh
```

### Mapping and quantification

* **Key libraries:** stringr

```
# map and quantify using Alevin-fry
scripts/Alevin-fry-unfiltered.sh
```

### merge metadata

* **Key libraries:** stringr

```
# merge XXX metadata
scripts/merge_metadata.sh
```

### QC

* **Key libraries:** stringr

```
# map and quantify using Alevin-fry
scripts/alevin-qc.sh
```

### Ambient RNA

* **Key libraries:** stringr

```
# run CellBender
scripts/cell-bender.sh
```

### Sample QC

* **Key libraries:** stringr

```
# remove doublets
scripts/mito-doublets.sh
```

### Integration, clustering, and filtering

* **Key libraries:** stringr

```
# iteratively integrate, cluster, and filter data
scripts/integration-clean.sh
```

### Annotation

* **Key libraries:** stringr

```
# prep allen_M1_ref
scripts/allen_M1_ref.sh
# prep allen_MCA_ref
allen_MCA_ref.sh
# prep allen_MTG_ref
allen_MTG_ref.sh
# prep BRAIN_ref
BRAIN_ref.sh
# prep BRAIN_ref_2
BRAIN_ref_2.sh
# prep data for annotation
prep-for-annotation.sh
# annotate
annotation-swarm.sh
# plot annotations
plot_azimuth.R
```

### pseudobulk samples

* **Key libraries:** stringr

```
# pseudobulk samples
scripts/pseudobulk.R
```

### test for sex diffs in subclass proportions

* **Key libraries:** stringr

```
# calculate and test proporitons
scripts/plot_calc_propotions.R
```

### differential expression

* **Key libraries:** stringr

```
# run differential expression (voom-limma)
scripts/differential_expression.sh
# run differential expression (mashr)
scripts/mashr.sh 
```

### transcriptome wide impact

* **Key libraries:** stringr

```
# estimate transcriptome wide impact for sex (sexTWI)
scripts/TRADE.sh
# leave-one-out
scripts/DE-loo-TRADE-jk.sh
# plot results
scripts/TRADE.R
```

### variance partitioning

* **Key libraries:** stringr

```
# run variance partitioning
scripts/variance-partitioning.sh  
# plot results
scripts/plot-variance-partitioning.R
```

### cluster autosomal sex effects

* **Key libraries:** stringr

```
# run homer
scripts/plot-sex-effects.R
```

### run homer TF enrichment on autosomal clusters

* **Key libraries:** stringr

```
# run homer
scripts/homer-clusters.sh  
# plot results
scripts/plot-variance-partitioning.R
```

### run GSEA on autosomal clusters

* **Key libraries:** stringr

```
# run GSEA
scripts/homer-clusters.sh  
# plot results
scripts/plot-variance-partitioning.R
```

### run sex-specific GWAS enrichments on autosomal clusters

* **Key libraries:** stringr

```
# download and prep GWAS data
scripts/homer-clusters.sh
# run enrichment
scripts/homer-clusters.sh  
# plot results
scripts/plot-variance-partitioning.R
```

### allele specific expression

* **Key libraries:** stringr

```
# estimate ASE
scripts/homer-clusters.sh  
# filter and plot results
scripts/plot-variance-partitioning.R
```


