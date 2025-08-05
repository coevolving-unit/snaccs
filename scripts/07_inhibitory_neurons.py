#!/usr/bin/env python3

"""
Inhibitory neuron subclustering and annotation
Specialized pipeline for analyzing inhibitory neurons
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
from subclustering_template import SubclusteringPipeline

def define_inhibitory_markers():
    """Define marker genes specific to inhibitory neuron subtypes"""
    marker_genes_dict = {
        'ExN': ['SLC17A7'],
        'CGE-derived': ['ADARB2'],
        'MGE-derived': ['LHX6'],
        'SST': ['SST'],
        'SNCG/PAX6': ['PAX6', 'SNCG'],
        'LAMP5': ['LAMP5'],
        'VIP': ['VIP'],
        'PVALB': ['PVALB'],
        'Chandelier': ['UNC5B'],
        'Reelin': ['RELN'],
        'Calretinin': ['CALB2'],
        'Cholecystokinin': ['CCK'],
        'NOS': ['NOS1', 'NOS2'],
        'NPY': ['NPY'],
        'ADAM33': ['ADAM33'],
        'Developmental': ['NDNF', 'ID2'],
    }
    return marker_genes_dict

def create_inhibitory_annotation_rounds():
    """Define cluster annotations for different rounds of analysis"""
    
    # Round 1: Quality filtering
    round1_annotation = {
        str(i): 'Keep' for i in range(33)
    }
    # Remove poor quality clusters
    round1_annotation.update({
        '29': 'Remove',  # high doublet
        '30': 'Remove',  # high doublet  
        '31': 'Remove',  # high umi high doublet
        '32': 'Remove',  # high umi high doublet
    })
    
    # Round 2: Initial cell type assignment
    round2_annotation = {
        '0': 'PVALB', '1': 'LAMP5', '2': 'VIP', '3': 'PVALB', '4': 'VIP',
        '5': 'PVALB', '6': 'SST', '7': 'VIP', '8': 'LAMP5', '9': 'SST',
        '10': 'VIP', '11': 'SST', '12': 'Chandelier', '13': 'LAMP5', '14': 'VIP',
        '15': 'PVALB', '16': 'SST', '17': 'SST', '18': 'VIP', '19': 'ADARB2',
        '20': 'SST', '21': 'PAX6', '22': 'ADARB2', '23': 'PAX6',
        '24': 'Remove',  # Mixed SST + ExN markers
    }
    
    # Round 3: Final refined annotation
    final_annotation = {
        '0': 'PVALB', '1': 'LAMP5', '2': 'PVALB', '3': 'VIP', '4': 'VIP',
        '5': 'SST', '6': 'SST', '7': 'LAMP5', '8': 'VIP', '9': 'VIP',
        '10': 'ADARB2', '11': 'SST', '12': 'Chandelier', '13': 'LAMP5', '14': 'VIP',
        '15': 'SST', '16': 'PVALB', '17': 'SST', '18': 'VIP', '19': 'PVALB',
        '20': 'PVALB', '21': 'SST', '22': 'PAX6', '23': 'VIP', '24': 'PAX6',
    }
    
    return round1_annotation, round2_annotation, final_annotation

class InhibitoryNeuronPipeline(SubclusteringPipeline):
    """Specialized pipeline for inhibitory neurons"""
    
    def __init__(self, input_file='adata_inn_corrected.h5ad'):
        super().__init__(input_file, 'InhibitoryNeurons', 'inn')
        
    def run_iterative_analysis(self, metadata_file, max_rounds=3):
        """Run iterative rounds of analysis with quality filtering"""
        
        print("Starting iterative inhibitory neuron analysis...")
        
        # Get annotation mappings
        round1_map, round2_map, final_map = create_inhibitory_annotation_rounds()
        
        # Round 1: Initial analysis and quality filtering
        print("\n=== ROUND 1: Initial Quality Assessment ===")
        self.load_data()
        self.process_data()
        self.integrate_data()
        self.cluster_data([1.0])
        self.add_metadata(metadata_file)
        
        # Create plots and find markers
        self.create_plots('leiden_res_1_0')
        self.find_marker_genes('leiden_res_1_0')
        
        # Annotate and filter
        self.annotate_clusters(round1_map, 'leiden_res_1_0')
        self.adata.write(f'{self.output_prefix}_round1.h5ad')
        
        # Filter to keep only good quality cells
        self.filter_cells_by_annotation(['Keep'])
        
        # Round 2: Initial cell type assignment
        print("\n=== ROUND 2: Initial Cell Type Assignment ===")
        self.reprocess_filtered_data()
        self.cluster_data([1.0])
        
        self.create_plots('leiden_res_1_0')
        self.find_marker_genes('leiden_res_1_0')
        
        # Create dotplot with markers
        markers = define_inhibitory_markers()
        self.create_dotplot(markers, 'leiden_res_1_0')
        
        # Annotate cell types
        self.annotate_clusters(round2_map, 'leiden_res_1_0')
        self.adata.write(f'{self.output_prefix}_round2.h5ad')
        
        # Filter to remove poor quality annotations
        keep_types = ['PVALB', 'LAMP5', 'VIP', 'SST', 'ADARB2', 'PAX6', 'Chandelier']
        self.filter_cells_by_annotation(keep_types)
        
        # Round 3: Final detailed annotation
        print("\n=== ROUND 3: Final Detailed Annotation ===")
        self.reprocess_filtered_data()
        self.cluster_data([0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
        
        self.create_plots('leiden_res_1_0')
        self.find_marker_genes('leiden_res_1_0')
        
        # Final annotation
        self.annotate_clusters(final_map, 'leiden_res_1_0')
        
        # Create final plots
        self.create_final_plots()
        
        # Save final results
        self.adata.write(f'{self.output_prefix}_final.h5ad')
        
        print("Inhibitory neuron analysis completed!")
        
    def filter_cells_by_annotation(self, keep_classes):
        """Filter cells based on annotation"""
        # Get raw data back
        raw_adata = self.adata.raw.to_adata()
        raw_adata.obs = self.adata.obs.copy()
        
        # Filter
        self.adata = raw_adata[raw_adata.obs['cell_class'].isin(keep_classes)]
        self.adata.raw = self.adata
        
        print(f"Filtered to {self.adata.shape[0]} cells")
        
    def reprocess_filtered_data(self):
        """Reprocess data after filtering"""
        print("Reprocessing filtered data...")
        self.process_data()
        self.integrate_data(max_iter=5)  # More harmony iterations for better integration
        
    def create_final_plots(self):
        """Create final publication-quality plots"""
        # Filter to final cell types
        final_types = ['PVALB', 'SST', 'VIP', 'PAX6', 'LAMP5', 'ADARB2', 'Chandelier']
        clean_adata = self.adata[self.adata.obs['cell_class'].isin(final_types)]
        
        # Create clean UMAP with custom colors
        sc.pl.umap(clean_adata, color=['cell_class'], 
                  palette=['#FF5733', '#FF4500', '#D95F02', '#FDCC8A', 
                          '#FFB347', '#FFC966', '#b45f06'],
                  legend_loc='on data', legend_fontsize=6, frameon=True,
                  save=f'umap_{self.output_prefix}_final_clean.png')
        
        # Create comprehensive dotplot
        marker_genes = {
            'CGE-derived': ['ADARB2'],
            'MGE-derived': ['LHX6'],
            'SST': ['SST'],
            'PAX6': ['PAX6'],
            'LAMP5': ['LAMP5'],
            'VIP': ['VIP'],
            'PVALB': ['PVALB'],
            'Chandelier': ['UNC5B'],
            'SNCG': ['SNCG'],
            'Reelin': ['RELN'],
            'Calretinin': ['CALB2'],
            'Cholecystokinin': ['CCK'],
            'NOS': ['NOS1', 'NOS2'],
            'NPY': ['NPY'],
            'ADAM33': ['ADAM33'],
            'Developmental': ['NDNF', 'ID2'],
        }
        
        sc.pl.dotplot(clean_adata, marker_genes, 'cell_class', 
                     dendrogram=True, save=f'dotplot_{self.output_prefix}_final.pdf')
        
        # Create additional marker expression UMAP
        key_markers = ['ADARB2', 'LHX6', 'SST', 'PAX6', 'SNCG', 'LAMP5', 'VIP', 'PVALB', 'UNC5B']
        sc.pl.umap(clean_adata, color=key_markers, legend_fontsize=6, frameon=True,
                  save=f'umap_{self.output_prefix}_markers.png')
        
        # Save cluster statistics
        cross_tab = pd.crosstab(clean_adata.obs['cell_class'], clean_adata.obs['orig.ident']).T
        cross_tab.to_csv(f'figures/{self.output_prefix}_final_sample_counts.csv')
        
        cross_tab_region = pd.crosstab(clean_adata.obs['cell_class'], clean_adata.obs['region']).T
        cross_tab_region.to_csv(f'figures/{self.output_prefix}_final_region_counts.csv')
        
    def create_supertype_annotation(self):
        """Create supertype annotations for finer clustering"""
        # This can be used for even more detailed subclustering
        supertype_mapping = {
            'PVALB1': 'PVALB', 'PVALB2': 'PVALB', 'PVALB3': 'PVALB',
            'PVALB4': 'PVALB', 'PVALB5': 'PVALB',
            'LAMP51': 'LAMP5', 'LAMP52': 'LAMP5', 'LAMP53': 'LAMP5',
            'VIP1': 'VIP', 'VIP2': 'VIP', 'VIP3': 'VIP', 'VIP4': 'VIP',
            'VIP5': 'VIP', 'VIP6': 'VIP', 'VIP7': 'VIP',
            'SST1': 'SST', 'SST2': 'SST', 'SST3': 'SST', 'SST4': 'SST',
            'SST5': 'SST', 'SST6': 'SST',
            'ADARB2': 'ADARB2',
            'PAX61': 'PAX6', 'PAX62': 'PAX6',
            'Chandelier': 'Chandelier',
        }
        return supertype_mapping


def main():
    """Main function for inhibitory neuron analysis"""
    # Initialize pipeline
    pipeline = InhibitoryNeuronPipeline()
    
    # Define metadata file
    metadata_file = '/data/region_samples_list.csv'
    
    # Run complete iterative analysis
    pipeline.run_iterative_analysis(metadata_file)
    
    print("Inhibitory neuron subclustering completed!")


if __name__ == "__main__":
    main()
