#!/usr/bin/env python3

"""
Single-cell RNA-seq data integration and initial clustering
Load and merge AnnData objects, perform integration with Harmony, and initial clustering
"""

import sys
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_sample_list(csv_file, exclude_bad_samples=True):
    """Load sample list and filter bad samples"""
    sample_list = pd.read_csv(csv_file).sample_id.tolist()
    
    # Bad samples identified from Cell Ranger QC and other issues
    bad_samples = [
        'sample_79536', 'sample_79544', 'sample_80953', 'sample_80956', 
        'sample_80957', 'sample_80958', 'sample_81687', 'sample_79531',
        'sample_79547', 'sample_79548', 'sample_80921', 'sample_80952', 
        'sample_80946'
    ]
    
    if exclude_bad_samples:
        sample_list = [x for x in sample_list if x not in bad_samples]
    
    # Remove 'sample_' prefix for filename matching
    sample_list = [sub.replace('sample_', '') for sub in sample_list]
    
    print(f"Found {len(sample_list)} valid samples")
    return sample_list

def merge_anndata_objects(sample_list):
    """Load and merge AnnData objects"""
    print("Loading AnnData objects...")
    
    filenames = [f"{s}_bender_filtered.h5ad" for s in sample_list]
    adatas = []
    
    for filename in filenames:
        try:
            adata = sc.read_h5ad(filename)
            adata.obs['orig.ident'] = filename
            adatas.append(adata)
            print(f"Loaded {filename}: {adata.shape}")
        except FileNotFoundError:
            print(f"Warning: {filename} not found, skipping")
    
    if not adatas:
        raise ValueError("No valid AnnData objects found")
    
    # Concatenate all objects
    adata_merged = adatas[0].concatenate(adatas[1:])
    print(f"Merged data shape: {adata_merged.shape}")
    
    return adata_merged

def add_metadata(adata, metadata_file):
    """Add region and sex information"""
    print("Adding metadata...")
    
    df = pd.read_csv(metadata_file, usecols=['sample_id', 'region', 'sex'])
    df['sample_id'] = df['sample_id'].str.replace("sample_", "")
    adata.obs['orig.ident'] = adata.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
    
    # Map region and sex
    region_dict = df.set_index('sample_id')['region'].to_dict()
    adata.obs['region'] = adata.obs['orig.ident'].map(region_dict).astype('category')
    
    sex_dict = df.set_index('sample_id')['sex'].to_dict()
    adata.obs['sex'] = adata.obs['orig.ident'].map(sex_dict).astype('category')
    
    return adata

def process_and_integrate(adata, max_iter_harmony=1):
    """Normalize, find variable genes, scale, PCA, and integrate with Harmony"""
    print("Processing data...")
    
    # Store raw data
    adata.raw = adata
    
    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]
    
    # Scale and PCA
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    
    print("Running Harmony integration...")
    # Harmony integration
    sce.pp.harmony_integrate(
        adata,
        key='orig.ident',
        basis='X_pca',
        adjusted_basis='X_pca_harmony',
        max_iter_harmony=max_iter_harmony
    )
    
    # Replace PCA with harmony-corrected PCA
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony'].copy()
    
    return adata

def cluster_and_umap(adata, resolution=1.0):
    """Perform clustering and UMAP"""
    print("Computing neighbors and UMAP...")
    
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, method='umap', metric='euclidean')
    
    # UMAP
    sc.tl.umap(adata, min_dist=0.5, spread=1.0, n_components=2)
    
    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution)
    
    return adata

def create_plots(adata, output_prefix=""):
    """Create standard QC and visualization plots"""
    print("Creating plots...")
    
    # UMAP plots
    sc.pl.umap(adata, color=['region'], palette=sns.color_palette("husl", 6),
              legend_fontsize=6, frameon=True, title='Region', 
              save=f'umap_regions_{output_prefix}.png')
    
    sc.pl.umap(adata, color=['sex'], palette=sns.color_palette("husl", 2),
              legend_fontsize=6, frameon=True, title='Sex',
              save=f'umap_sex_{output_prefix}.png')
    
    sc.pl.umap(adata, color=['orig.ident'], palette=sns.color_palette("husl", 200),
              legend_fontsize=6, frameon=True, title='Donor',
              save=f'umap_donors_{output_prefix}.png')
    
    sc.pl.umap(adata, color=['leiden'], palette=sns.color_palette("husl", 36),
              legend_loc='on data', legend_fontsize=6, frameon=True, title='Leiden',
              save=f'umap_leiden_{output_prefix}.png')
    
    # QC metrics
    sc.pl.umap(adata, color=['percent.mt'], legend_fontsize=6, frameon=True, 
              title='Mitochondrial %', save=f'umap_mito_{output_prefix}.png')
    
    sc.pl.umap(adata, color=['nCount_RNA'], legend_fontsize=6, frameon=True, 
              title='UMI Count', vmax=30000, save=f'umap_umi_{output_prefix}.png')
    
    # Violin plots
    sc.pl.violin(adata, keys='percent.mt', groupby='leiden', rotation=90,
                save=f'violin_mito_{output_prefix}.png')
    
    sc.pl.violin(adata, keys='nCount_RNA', groupby='leiden', rotation=90,
                save=f'violin_umi_{output_prefix}.png')

def save_cluster_stats(adata, output_prefix=""):
    """Save cluster composition statistics"""
    print("Saving cluster statistics...")
    
    # Sample counts per cluster
    cross_tab = pd.crosstab(adata.obs['leiden'], adata.obs['orig.ident']).T
    cross_tab.to_csv(f'figures/{output_prefix}_sample_counts.csv')
    
    # Sample proportions per cluster
    cross_tab_prop = pd.crosstab(adata.obs['leiden'], adata.obs['orig.ident'], 
                                normalize='columns').T
    cross_tab_prop.to_csv(f'figures/{output_prefix}_sample_proportions.csv')

def main():
    """Main processing pipeline"""
    print("Starting single-cell RNA-seq integration pipeline...")
    
    # Load sample list
    sample_list = load_sample_list('region_samples_list.csv')
    
    # Merge data
    adata = merge_anndata_objects(sample_list)
    
    # Add metadata
    adata = add_metadata(adata, '/data/DNU/alex/snRNAseq/alevin-all/region_samples_list.csv')
    
    # Save raw merged data
    adata.write('adata_corrected.h5ad')
    
    # Process and integrate
    adata = process_and_integrate(adata, max_iter_harmony=1)
    adata.write('adata_harmony_corrected.h5ad')
    
    # Clustering and UMAP
    adata = cluster_and_umap(adata, resolution=1.0)
    adata.write('adata_leiden_res1_corrected.h5ad')
    
    # Create plots
    create_plots(adata, "corrected")
    
    # Save statistics
    save_cluster_stats(adata, "R1_corrected")
    
    print("Integration pipeline completed successfully!")

if __name__ == "__main__":
    main()
