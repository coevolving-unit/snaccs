#!/usr/bin/env python3

"""
Cell type annotation and marker gene analysis
Annotate clusters based on marker genes and create dotplots
"""

import scanpy as sc
import pandas as pd
import seaborn as sns

# Set scanpy settings
sc.settings.verbosity = 3

def define_marker_genes():
    """Define marker genes for different cell types"""
    marker_genes_dict = {
        'Oligo': ['PLP1', 'OLIG1', 'OLIG2', 'MOBP', 'BCAS1', 'OPALIN'],
        'OPC': ['SOX10', 'CSPG4', 'PDGFRA', 'LUZP2', 'PCDH15'],
        'Micro': ['CD74', 'APBB1IP', 'COBLL1', 'EBF1'],
        'Macrophages': ['MRC1'],
        'T cells': ['PTPRC'],
        'Astro': ['GFAP', 'ALDH1L1', 'AQP4', 'GJA1'],
        'Neurons': ['GRIN1', 'GRIN2B'],
        'ExN': ['NEUROD6', 'NRGN', 'SLC17A6', 'SLC17A7', 'SATB2'],
        'InN': ['GAD1', 'GAD2', 'CALB2', 'PVALB', 'VIP', 'SST', 'LAMP5'],
        'Pericytes': ['MCAM'],
        'Endo': ['FLT1'],
        'Fibro': ['DCN'],
        'Ambient RNA': ['SYT1', 'CSMD1', 'KCNIP4', 'RBFOX1', 'RALYL', 
                       'NRGN', 'CHN1', 'MT-ND1', 'MT-CO1'],
    }
    return marker_genes_dict

def create_cluster_annotation_mapping():
    """Define cluster to cell type annotation mapping"""
    # This mapping should be updated based on manual inspection of markers
    cluster2annotation = {
        '0': 'ExN',    # Excitatory neurons
        '1': 'Non-neuronal',  # Oligodendrocytes
        '2': 'ExN',
        '3': 'ExN',
        '4': 'ExN',
        '5': 'ExN',
        '6': 'InN',    # Inhibitory neurons
        '7': 'InN',
        '8': 'InN',
        '9': 'Remove',  # Low quality - high mito
        '10': 'ExN',
        '11': 'ExN',
        '12': 'Non-neuronal',  # OPC
        '13': 'ExN',
        '14': 'Non-neuronal',  # Astrocytes
        '15': 'InN',
        '16': 'ExN',
        '17': 'ExN',
        '18': 'Non-neuronal',  # Microglia
        '19': 'InN',
        '20': 'ExN',
        '21': 'Non-neuronal',  # Pericytes
        '22': 'ExN',
        '23': 'ExN',
        '24': 'ExN',
        '25': 'Non-neuronal',  # Pericytes
        '26': 'InN',
        '27': 'InN',
        '28': 'InN',
        '29': 'Non-neuronal',  # Pericytes
        '30': 'ExN',
        '31': 'Non-neuronal',  # T cells
        '32': 'Remove',  # Mixed markers
        '33': 'ExN',
        '34': 'Remove',  # No clear markers
        '35': 'Remove',  # Strange markers
        '36': 'Remove',  # No markers
        '37': 'Remove',
        '38': 'Remove',
        '39': 'Remove',
        '40': 'Remove',
    }
    return cluster2annotation

def annotate_clusters(adata, cluster_col='leiden', annotation_mapping=None):
    """Annotate clusters with cell type labels"""
    if annotation_mapping is None:
        annotation_mapping = create_cluster_annotation_mapping()
    
    adata.obs['cell_class'] = adata.obs[cluster_col].map(annotation_mapping).astype('category')
    
    print("Cell class counts:")
    print(adata.obs['cell_class'].value_counts())
    
    return adata

def create_dotplot(adata, marker_genes_dict, groupby='leiden', save_name='dotplot.png'):
    """Create dotplot of marker genes"""
    sc.pl.dotplot(adata, marker_genes_dict, groupby, 
                 dendrogram=True, save=save_name)

def subset_by_cell_class(adata, cell_classes, output_file):
    """Subset data by cell class and save"""
    # Get original data with all genes
    original_adata = adata.raw.to_adata()
    original_adata.obs = adata.obs.copy()
    
    # Subset by cell classes
    subset_adata = original_adata[original_adata.obs['cell_class'].isin(cell_classes)]
    
    print(f"Subset {cell_classes} shape: {subset_adata.shape}")
    subset_adata.write(output_file)
    
    return subset_adata

def main():
    """Main annotation pipeline"""
    print("Starting cell type annotation...")
    
    # Load integrated data
    adata = sc.read_h5ad('adata_leiden_res1_corrected.h5ad')
    
    # Define marker genes
    marker_genes = define_marker_genes()
    
    # Create dotplot before annotation
    create_dotplot(adata, marker_genes, 'leiden', 'R1_dotplot_corrected.png')
    
    # Annotate clusters
    adata = annotate_clusters(adata)
    
    # Plot cell classes
    sc.pl.umap(adata, color=['cell_class'], palette=sns.color_palette("husl", 4),
              legend_fontsize=6, legend_loc='on data', frameon=True, 
              title='Cell Class', save='cell_class_corrected.png')
    
    # Save annotated data
    adata.write('adata_cell_classes_corrected.h5ad')
    
    # Create subsets for different cell types
    print("Creating cell type subsets...")
    
    # Excitatory neurons
    subset_by_cell_class(adata, ['ExN'], 'adata_exn_corrected.h5ad')
    
    # Inhibitory neurons  
    subset_by_cell_class(adata, ['InN'], 'adata_inn_corrected.h5ad')
    
    # Non-neuronal cells
    subset_by_cell_class(adata, ['Non-neuronal'], 'adata_nonneuronal_corrected.h5ad')
    
    print("Cell type annotation completed!")

if __name__ == "__main__":
    main()
