#!/usr/bin/env python3

"""
Final integration and visualization
Combine all cell type analyses and create final integrated dataset
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def merge_cell_type_annotations():
    """Merge annotations from all cell type analyses"""
    print("Loading individual cell type analyses...")
    
    # Load each analyzed cell type
    other = sc.read_h5ad('other_subtype_final.h5ad')
    inn = sc.read_h5ad('inn_subtype_final.h5ad')  
    exn = sc.read_h5ad('exn_subtype_final.h5ad')
    
    print(f"Non-neuronal cells: {other.shape}")
    print(f"Inhibitory neurons: {inn.shape}")  
    print(f"Excitatory neurons: {exn.shape}")
    
    # Extract cell annotations
    other_df = pd.DataFrame({
        'cell_class': other.obs['cell_class'],
        'supertype': other.obs.get('supertype', other.obs['cell_class']),
        'cell_id': other.obs.index
    })
    
    inn_df = pd.DataFrame({
        'cell_class': inn.obs['cell_class'],
        'supertype': inn.obs.get('supertype', inn.obs['cell_class']),
        'cell_id': inn.obs.index
    })
    
    exn_df = pd.DataFrame({
        'cell_class': exn.obs['cell_class'],
        'supertype': exn.obs.get('supertype', exn.obs['cell_class']),
        'cell_id': exn.obs.index
    })
    
    # Combine annotations
    all_annotations = pd.concat([other_df, inn_df, exn_df], ignore_index=True)
    all_annotations = all_annotations.set_index('cell_id')
    
    print(f"Total annotated cells: {len(all_annotations)}")
    
    return all_annotations

def create_integrated_dataset(annotations_df):
    """Create final integrated dataset with all annotations"""
    print("Creating integrated dataset...")
    
    # Load original integrated data
    adata_orig = sc.read_h5ad('adata_corrected.h5ad')
    
    # Filter to cells with annotations
    cell_ids = annotations_df.index.tolist()
    adata_integrated = adata_orig[adata_orig.obs.index.isin(cell_ids)]
    
    # Add annotations
    adata_integrated.obs['cell_class'] = annotations_df.loc[adata_integrated.obs.index, 'cell_class']
    adata_integrated.obs['supertype'] = annotations_df.loc[adata_integrated.obs.index, 'supertype']
    
    # Add metadata if not present
    if 'region' not in adata_integrated.obs.columns:
        add_metadata(adata_integrated)
    
    print(f"Integrated dataset shape: {adata_integrated.shape}")
    print("\nCell type counts:")
    print(adata_integrated.obs['cell_class'].value_counts())
    
    return adata_integrated

def add_metadata(adata):
    """Add sample metadata to integrated dataset"""
    print("Adding metadata...")
    
    metadata_file = '/data/DNU/alex/snRNAseq/alevin-all/region_samples_list.csv'
    df = pd.read_csv(metadata_file, usecols=['sample_id', 'region', 'sex', 'batch', 'new_id'])
    df['sample_id'] = df['sample_id'].str.replace("sample_", "")
    df['new_id'] = df['new_id'].str.replace("sample_", "")
    
    # Clean orig.ident column
    adata.obs['orig.ident'] = adata.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
    
    # Map metadata
    for col in ['region', 'sex', 'batch', 'new_id']:
        col_dict = df.set_index('sample_id')[col].to_dict()
        adata.obs[col] = adata.obs['orig.ident'].map(col_dict).astype('category')

def recompute_integration(adata):
    """Recompute integration for final dataset"""
    print("Recomputing integration for final dataset...")
    
    # Store raw data
    adata.raw = adata
    
    # Standard processing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Harmony integration
    import scanpy.external as sce
    sce.pp.harmony_integrate(
        adata,
        key='orig.ident',
        basis='X_pca',
        adjusted_basis='X_pca_harmony',
        max_iter_harmony=3
    )
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony'].copy()
    
    # Compute neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.5, spread=1.0)
    
    return adata

def create_final_visualizations(adata):
    """Create comprehensive final visualizations"""
    print("Creating final visualizations...")
    
    # Define color palettes for different cell types
    cell_type_colors = {
        # Excitatory neurons - greens
        'L23IT': '#66C2A5', 'L4IT': '#0C7C59', 'L5IT': '#99C749', 
        'L6IT': '#78C679', 'L5ET': '#C2E699', 'L6CT': '#4EB265', 
        'L6b': '#a2c4c9', 'L56NP': '#1B9E77',
        
        # Inhibitory neurons - oranges/reds  
        'PVALB': '#FF5733', 'SST': '#FF4500', 'VIP': '#D95F02',
        'PAX6': '#FDCC8A', 'LAMP5': '#FFB347', 'ADARB2': '#FFC966',
        'Chandelier': '#b45f06',
        
        # Non-neuronal - purples/blues
        'Oligo': '#800080', 'OPC': '#9370DB', 'Astro': '#6A5ACD',
        'Micro': '#D8BFD8', 'Tcell': '#7E78B8', 'BEC': '#8E89A0',
        'Fibro': '#BCBDDC', 'Peri': '#8c5c86', 'SM': '#a64d79'
    }
    
    # Main UMAP plots
    sc.pl.umap(adata, color=['cell_class'], 
              palette=[cell_type_colors.get(ct, '#888888') for ct in adata.obs['cell_class'].cat.categories],
              legend_loc='on data', legend_fontsize=6, frameon=True,
              title='Cell Types', save='final_integrated_celltypes.png')
    
    sc.pl.umap(adata, color=['region'], palette=sns.color_palette("husl", 6),
              legend_fontsize=6, frameon=True, title='Brain Region',
              save='final_integrated_regions.png')
    
    sc.pl.umap(adata, color=['sex'], palette=sns.color_palette("husl", 2),
              legend_fontsize=6, frameon=True, title='Sex',
              save='final_integrated_sex.png')
    
    # QC metrics
    sc.pl.umap(adata, color=['percent.mt'], legend_fontsize=6, frameon=True,
              title='Mitochondrial %', save='final_integrated_mito.png')
    
    sc.pl.umap(adata, color=['nCount_RNA'], legend_fontsize=6, frameon=True,
              title='UMI Count', vmax=30000, save='final_integrated_umi.png')
    
    # Create summary dotplot with key markers
    create_summary_dotplot(adata)
    
    # Create composition plots
    create_composition_plots(adata)

def create_summary_dotplot(adata):
    """Create comprehensive dotplot with all cell type markers"""
    
    all_markers = {
        # Excitatory neurons
        'ExN': ['SLC17A7', 'SATB2'],
        'IT': ['RORB'],
        'Upper IT': ['CUX2'],
        'Lower IT': ['THEMIS'],
        'L5 ET': ['ADRA1A'],
        'L6 CT': ['MCC'],
        'L6b': ['DPP4'],
        'L5/6 NP': ['ZNF385D'],
        
        # Inhibitory neurons
        'InN': ['GAD1', 'GAD2'],
        'PVALB': ['PVALB'],
        'SST': ['SST'],
        'VIP': ['VIP'],
        'LAMP5': ['LAMP5'],
        'PAX6': ['PAX6'],
        'CGE': ['ADARB2'],
        'Chandelier': ['UNC5B'],
        
        # Non-neuronal
        'Astro': ['AQP4', 'GFAP'],
        'Oligo': ['MOG', 'PLP1'],
        'OPC': ['PDGFRA', 'SOX10'],
        'Micro': ['APBB1IP'],
        'Endo': ['FLT1'],
        'Peri': ['GRM8'],
        'T cells': ['PTPRC'],
    }
    
    sc.pl.dotplot(adata, all_markers, 'cell_class', dendrogram=True,
                 save='final_integrated_markers.pdf')

def create_composition_plots(adata):
    """Create cell type composition analysis plots"""
    
    # Sample composition
    sample_comp = pd.crosstab(adata.obs['cell_class'], adata.obs['orig.ident'], normalize='columns')
    sample_comp.to_csv('figures/final_sample_composition.csv')
    
    # Region composition  
    region_comp = pd.crosstab(adata.obs['cell_class'], adata.obs['region'], normalize='columns')
    region_comp.to_csv('figures/final_region_composition.csv')
    
    # Sex composition
    sex_comp = pd.crosstab(adata.obs['cell_class'], adata.obs['sex'], normalize='columns')
    sex_comp.to_csv('figures/final_sex_composition.csv')
    
    # Create stacked bar plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Region composition plot
    region_comp.T.plot(kind='bar', stacked=True, ax=axes[0], 
                      colormap='tab20', legend=False)
    axes[0].set_title('Cell Type Composition by Region')
    axes[0].set_xlabel('Brain Region')
    axes[0].set_ylabel('Proportion')
    
    # Sex composition plot  
    sex_comp.T.plot(kind='bar', stacked=True, ax=axes[1],
                   colormap='tab20', legend=False)
    axes[1].set_title('Cell Type Composition by Sex')
    axes[1].set_xlabel('Sex')
    axes[1].set_ylabel('Proportion')
    
    # Overall cell type counts
    cell_counts = adata.obs['cell_class'].value_counts()
    cell_counts.plot(kind='bar', ax=axes[2], color='skyblue')
    axes[2].set_title('Total Cell Counts by Type')
    axes[2].set_xlabel('Cell Type')
    axes[2].set_ylabel('Number of Cells')
    
    plt.tight_layout()
    plt.savefig('figures/final_composition_summary.png', dpi=300
