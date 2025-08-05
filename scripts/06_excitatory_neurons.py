#!/usr/bin/env python3

"""
Excitatory neuron subclustering and annotation
Specialized pipeline for analyzing excitatory neurons
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
from subclustering_template import SubclusteringPipeline

def define_excitatory_markers():
    """Define marker genes specific to excitatory neuron subtypes"""
    marker_genes_dict = {
        'Neurons': ['GRIN1', 'GRIN2B'],
        'ExN': ['NEUROD6', 'NRGN', 'SLC17A6', 'SLC17A7', 'SATB2'],
        'IT': ['RORB'],
        'upper-layer IT': ['CUX2', 'LAMP5'],
        'lower-layer IT': ['THEMIS'],
        'lower-layer non-IT': ['FEZF2'],
        'L5 ET': ['ADRA1A'],
        'L6 CT': ['MCC'],
        'L6b': ['DPP4'],
        'L5/6 NP': ['ZNF385D'],
        'Ambient RNA': ['SYT1', 'CSMD1', 'KCNIP4', 'RALYL', 'NRGN', 'CHN1', 'MT-ND1', 'MT-CO1'],
    }
    return marker_genes_dict

def create_excitatory_annotation_rounds():
    """Define cluster annotations for different rounds of analysis"""
    
    # Round 1: Initial quality filtering
    round1_annotation = {
        '0': 'Keep', '1': 'Keep', '2': 'Keep', '3': 'Keep', '4': 'Keep',
        '5': 'Keep', '6': 'Keep', '7': 'Keep', '8': 'Keep', '9': 'Keep',
        '10': 'Keep', '11': 'Keep', '12': 'Keep', '13': 'Keep', '14': 'Keep',
        '15': 'Keep', '16': 'Keep', '17': 'Keep', '18': 'Keep',
        '19': 'Remove',  # high mito + high doublet + low intronic + high UMI
        '20': 'Remove',  # donor bias + high mito
        '21': 'Keep',    # low markers + low UMI (all samples have)
        '22': 'Remove',  # donor bias + high doublet
        '23': 'Keep', '24': 'Keep',
        '25': 'Remove',  # mixed markers w/ oligo + high doublet
    }
    
    # Round 2: Layer assignment
    round2_annotation = {
        '0': 'upper', '1': 'lower', '2': 'upper', '3': 'upper', '4': 'upper',
        '5': 'upper', '6': 'lower', '7': 'lower', '8': 'upper', '9': 'upper',
        '10': 'upper', '11': 'upper', '12': 'upper', '13': 'lower', '14': 'upper',
        '15': 'lower', '16': 'lower', '17': 'upper', '18': 'lower', '19': 'lower',
        '20': 'upper', '21': 'lower', '22': 'lower',
        '23': 'Remove',  # low markers low UMI low intronic
        '24': 'Remove', '25': 'Remove', '26': 'Remove',
        '27': 'lower', '28': 'upper', '29': 'lower', '30': 'upper',
    }
    
    # Round 4: Final detailed annotation
    final_annotation = {
        '0': 'L6IT', '1': 'L23IT', '2': 'L5IT', '3': 'L23IT', '4': 'L23IT',
        '5': 'L23IT', '6': 'L23IT', '7': 'L23IT', '8': 'L4IT', '9': 'L5IT',
        '10': 'L4IT', '11': 'L23IT', '12': 'L23IT', '13': 'L6CT', '14': 'L23IT',
        '15': 'L5IT', '16': 'L23IT', '17': 'L23IT', '18': 'L4IT', '19': 'L6IT',
        '20': 'L5IT', '21': 'L4IT', '22': 'L56NP', '23': 'L5IT', '24': 'L6b',
        '25': 'L23IT', '26': 'L6b', '27': 'L5ET', '28': 'Remove', '29': 'Remove',
        '30': 'Remove', '31': 'Remove', '32': 'L23IT',
    }
    
    return round1_annotation, round2_annotation, final_annotation

class ExcitatoryNeuronPipeline(SubclusteringPipeline):
    """Specialized pipeline for excitatory neurons"""
    
    def __init__(self, input_file='adata_exn_corrected.h5ad'):
        super().__init__(input_file, 'ExcitatoryNeurons', 'exn')
        
    def run_iterative_analysis(self, metadata_file, max_rounds=4):
        """Run iterative rounds of analysis with quality filtering"""
        
        print("Starting iterative excitatory neuron analysis...")
        
        # Get annotation mappings
        round1_map, round2_map, final_map = create_excitatory_annotation_rounds()
        
        # Round 1: Initial analysis and quality filtering
        print("\n=== ROUND 1: Initial Quality Assessment ===")
        self.load_data()
        self.process_data()
        self.integrate_data()
        self.cluster_data([1.8])  # Use resolution 1.8 as in original
        self.add_metadata(metadata_file)
        
        # Create plots and find markers
        self.create_plots('leiden_res_1_8')
        self.find_marker_genes('leiden_res_1_8')
        
        # Annotate and filter
        self.annotate_clusters(round1_map, 'leiden_res_1_8')
        self.adata.write(f'{self.output_prefix}_round1.h5ad')
        
        # Filter to keep only good quality cells
        self.filter_cells_by_annotation(['Keep'])
        
        # Round 2: Layer assignment
        print("\n=== ROUND 2: Layer Assignment ===")
        self.reprocess_filtered_data()
        self.cluster_data([1.8])
        
        self.create_plots('leiden_res_1_8')
        self.find_marker_genes('leiden_res_1_8')
        
        # Annotate layers
        self.annotate_clusters(round2_map, 'leiden_res_1_8')
        self.adata.write(f'{self.output_prefix}_round2.h5ad')
        
        # Filter again
        self.filter_cells_by_annotation(['upper', 'lower'])
        
        # Round 3: Quality check
        print("\n=== ROUND 3: Final Quality Check ===")
        self.reprocess_filtered_data()
        self.cluster_data([1.8])
        
        # Simple keep/remove annotation (implement based on QC metrics)
        round3_map = {str(i): 'keep' for i in range(40)}
        round3_map.update({'14': 'Remove', '28': 'Remove'})  # Low quality clusters
        
        self.annotate_clusters(round3_map, 'leiden_res_1_8')
        self.filter_cells_by_annotation(['keep'])
        
        # Round 4: Final detailed annotation
        print("\n=== ROUND 4: Final Detailed Annotation ===")
        self.reprocess_filtered_data()
        self.cluster_data([1.0, 1.4, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 2.9, 3.0])
        
        self.create_plots('leiden_res_1_8')
        self.find_marker_genes('leiden_res_1_8')
        
        # Final annotation
        self.annotate_clusters(final_map, 'leiden_res_1_8')
        
        # Create final plots with clean cell types
        self.create_final_plots()
        
        # Save final results
        self.adata.write(f'{self.output_prefix}_final.h5ad')
        
        print("Excitatory neuron analysis completed!")
        
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
        self.integrate_data()
        
    def create_final_plots(self):
        """Create final publication-quality plots"""
        # Filter to final cell types
        final_types = ['L23IT', 'L4IT', 'L5IT', 'L6IT', 'L5ET', 'L6CT', 'L6b', 'L56NP']
        clean_adata = self.adata[self.adata.obs['cell_class'].isin(final_types)]
        
        # Create clean UMAP
        sc.pl.umap(clean_adata, color=['cell_class'], 
                  palette=['#66C2A5', '#0C7C59', '#99C749', '#78C679', 
                          '#C2E699', '#4EB265', '#a2c4c9', '#1B9E77'],
                  legend_loc='on data', legend_fontsize=6, frameon=True,
                  save=f'umap_{self.output_prefix}_final_clean.png')
        
        # Create dotplot with key markers
        marker_genes = define_excitatory_markers()
        sc.pl.dotplot(clean_adata, marker_genes, 'cell_class', 
                     dendrogram=True, save=f'dotplot_{self.output_prefix}_final.pdf')
        
        # Save cluster statistics
        cross_tab = pd.crosstab(clean_adata.obs['cell_class'], clean_adata.obs['orig.ident']).T
        cross_tab.to_csv(f'figures/{self.output_prefix}_final_sample_counts.csv')
        
        cross_tab_region = pd.crosstab(clean_adata.obs['cell_class'], clean_adata.obs['region']).T
        cross_tab_region.to_csv(f'figures/{self.output_prefix}_final_region_counts.csv')


def main():
    """Main function for excitatory neuron analysis"""
    # Initialize pipeline
    pipeline = ExcitatoryNeuronPipeline()
    
    # Define metadata file
    metadata_file = '/data/region_samples_list.csv'
    
    # Run complete iterative analysis
    pipeline.run_iterative_analysis(metadata_file)
    
    print("Excitatory neuron subclustering completed!")


if __name__ == "__main__":
    main()
