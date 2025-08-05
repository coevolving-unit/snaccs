#!/usr/bin/env python3

"""
Template for subclustering analysis
This script provides a template for detailed subclustering of specific cell types
Modify the cell_type_markers and cluster_annotations for specific cell types
"""

import sys
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sns

# Set scanpy settings
sc.settings.verbosity = 3

class SubclusteringPipeline:
    def __init__(self, input_file, cell_type, output_prefix):
        self.input_file = input_file
        self.cell_type = cell_type
        self.output_prefix = output_prefix
        self.adata = None
        
    def load_data(self):
        """Load the subset data"""
        print(f"Loading {self.cell_type} data from {self.input_file}")
        self.adata = sc.read_h5ad(self.input_file)
        self.adata.raw = self.adata
        print(f"Data shape: {self.adata.shape}")
        
    def process_data(self, target_sum=1e4, n_top_genes=2000):
        """Standard preprocessing pipeline"""
        print("Processing data...")
        
        # Normalize and log transform
        sc.pp.normalize_total(self.adata, target_sum=target_sum)
        sc.pp.log1p(self.adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(self.adata, n_top_genes=n_top_genes)
        self.adata = self.adata[:, self.adata.var.highly_variable]
        
        # Scale data
        sc.pp.scale(self.adata, max_value=10)
        
        # PCA
        sc.tl.pca(self.adata, svd_solver='arpack')
        
    def integrate_data(self, max_iter=1):
        """Harmony integration"""
        print("Running Harmony integration...")
        sce.pp.harmony_integrate(
            self.adata,
            key='orig.ident',
            basis='X_pca',
            adjusted_basis='X_pca_harmony',
            max_iter_harmony=max_iter
        )
        self.adata.obsm['X_pca'] = self.adata.obsm['X_pca_harmony'].copy()
        
    def cluster_data(self, resolutions=[0.4, 0.6, 0.8, 1.0, 1.2, 1.4]):
        """Clustering at multiple resolutions"""
        print("Computing neighbors and clustering...")
        
        # Compute neighborhood graph
        sc.pp.neighbors(self.adata, n_neighbors=15, n_pcs=50, 
                       method='umap', metric='euclidean')
        
        # UMAP
        sc.tl.umap(self.adata, min_dist=0.5, spread=1.0, n_components=2)
        
        # Leiden clustering at multiple resolutions
        for res in resolutions:
            res_key = f"leiden_res_{str(res).replace('.', '_')}"
            sc.tl.leiden(self.adata, resolution=res, key_added=res_key)
            print(f"Resolution {res}: {len(self.adata.obs[res_key].cat.categories)} clusters")
    
    def add_metadata(self, metadata_file):
        """Add sample metadata"""
        print("Adding metadata...")
        
        df = pd.read_csv(metadata_file, usecols=['sample_id', 'region', 'sex', 'batch'])
        df['sample_id'] = df['sample_id'].str.replace("sample_", "")
        self.adata.obs['orig.ident'] = self.adata.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
        
        # Map metadata
        for col in ['region', 'sex', 'batch']:
            col_dict = df.set_index('sample_id')[col].to_dict()
            self.adata.obs[col] = self.adata.obs['orig.ident'].map(col_dict).astype('category')
    
    def find_marker_genes(self, groupby='leiden_res_1_0', n_genes=20, method='wilcoxon'):
        """Find marker genes for clusters"""
        print("Finding marker genes...")
        
        # Subsample for speed if dataset is large
        if self.adata.n_obs > 20000:
            target_cells = 200
            subset_data = []
            
            for cluster in self.adata.obs[groupby].cat.categories:
                cluster_data = self.adata[self.adata.obs[groupby] == cluster]
                if cluster_data.n_obs > target_cells:
                    sc.pp.subsample(cluster_data, n_obs=target_cells)
                subset_data.append(cluster_data)
            
            adata_sub = subset_data[0].concatenate(subset_data[1:])
        else:
            adata_sub = self.adata
        
        # Find markers
        sc.tl.rank_genes_groups(adata_sub, groupby, method=method, n_genes=n_genes)
        
        # Create results dataframe
        result = adata_sub.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        markers_df = pd.DataFrame({
            group + '_' + key[:1]: result[key][group]
            for group in groups 
            for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']
        })
        
        # Save results
        markers_df.to_csv(f'{self.output_prefix}_marker_genes.csv')
        
        # Plot top markers
        sc.pl.rank_genes_groups(adata_sub, n_genes=10, sharey=False,
                               save=f'{self.output_prefix}_marker_genes.png')
    
    def create_plots(self, cluster_col='leiden_res_1_0'):
        """Create standard visualization plots"""
        print("Creating plots...")
        
        plot_params = [
            ('region', sns.color_palette("husl", 6), 'Region'),
            ('sex', sns.color_palette("husl", 2), 'Sex'),
            ('orig.ident', sns.color_palette("husl", 200), 'Donor'),
            (cluster_col, sns.color_palette("husl", 36), 'Clusters')
        ]
        
        for col, palette, title in plot_params:
            sc.pl.umap(self.adata, color=[col], palette=palette,
                      legend_fontsize=6, frameon=True, title=title,
                      save=f'umap_{self.output_prefix}_{col}.png')
        
        # QC metrics
        sc.pl.umap(self.adata, color=['percent.mt'], legend_fontsize=6, 
                  frameon=True, title='Mitochondrial %',
                  save=f'umap_{self.output_prefix}_mito.png')
        
        sc.pl.umap(self.adata, color=['nCount_RNA'], legend_fontsize=6,
                  frameon=True, title='UMI Count', vmax=30000,
                  save=f'umap_{self.output_prefix}_umi.png')
        
        # Violin plots
        sc.pl.violin(self.adata, keys='percent.mt', groupby=cluster_col, 
                    rotation=90, save=f'violin_{self.output_prefix}_mito.png')
        
        sc.pl.violin(self.adata, keys='nCount_RNA', groupby=cluster_col,
                    rotation=90, save=f'violin_{self.output_prefix}_umi.png')
    
    def save_cluster_stats(self, cluster_col='leiden_res_1_0'):
        """Save cluster statistics"""
        print("Saving cluster statistics...")
        
        # Sample counts and proportions
        cross_tab = pd.crosstab(self.adata.obs[cluster_col], self.adata.
