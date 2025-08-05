#!/usr/bin/env python3

"""
Non-neuronal cell subclustering and annotation
Specialized pipeline for analyzing glial and vascular cells
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
from subclustering_template import SubclusteringPipeline

def define_nonneuronal_markers():
    """Define marker genes specific to non-neuronal cell subtypes"""
    marker_genes_dict = {
        'Astro': ['AQP4', 'GFAP'],
        'OPC': ['PDGFRA', 'SOX10', 'BCAN', 'SOX6'],
        'COP': ['EPHB1', 'SH3RF3', 'GPR17', 'SOX6'],
        'Oligo': ['MOG', 'PLP1', 'MOBP', 'OPALIN', 'KLK6', 'OLIG1', 'OLIG2', 'ITPR2'],
        'Immune oligo': ['APOE', 'CD74'],
        'Endo': ['CLDN5', 'FLT1'],
        'Fibro': ['DCN'],
        'Peri': ['MUSTN1', 'MCAM'],
        'BEC-cap': ['FLT1'],
        'BEC-venous': ['TSHZ2'],
        'BEC-arterial': ['ARL15'],
        'Ependymal': ['CFAP299'],
        'Meningeal Fibroblast': ['SLC4A4'],
        'Perivascular Fibroblast': ['CEMIP'],
        'Pericyte': ['GRM8'],
        'Smooth muscle': ['SLIT3'],
        'Micro': ['APBB1IP'],
        'Myeloid microglia': ['AIF1', 'CD14'],
        'Macrophage': ['MRC1'],
        'Tcells': ['PTPRC'],
    }
    return marker_genes_dict

def create_nonneuronal_annotation_rounds():
    """Define cluster annotations for different rounds of analysis"""
    
    # Round 1: Initial broad classification
    round1_annotation = {
        '0': 'Oligo', '1': 'Oligo', '2': 'Astro', '3': 'OPC', '4': 'OPC',
        '5': 'Oligo', '6': 'Oligo', '7': 'Immune', '8': 'Endo', '9': 'Astro',
        '10': 'Endo', '11': 'Oligo', '12': 'Endo', '13': 'Endo', '14': 'Immune',
        '15': 'Immune', '16': 'Endo', '17': 'Remove', '18': 'Remove', '19': 'Remove',
        '20': 'Remove', '21': 'Astro', '22': 'Endo', '23': 'Oligo', '24': 'Endo',
        '25': 'Remove', '26': 'OPC', '27': 'Remove', '28': 'Oligo',
    }
    
    # Round 2: More detailed classification
    round2_annotation = {
        '0': 'Astro', '1': 'Oligo', '2': 'Oligo', '3': 'Oligo', '4': 'OPC',
        '5': 'Oligo', '6': 'OPC', '7': 'Micro', '8': 'BEC-cap', '9': 'Remove',
        '10': 'Astro', '11': 'BEC-arterial', '12': 'Perivasc Fibro', '13': 'Remove',
        '14': 'OPC', '15': 'Oligo', '16': 'T cells', '17': 'Pericyte', '18': 'BEC-venous',
        '19': 'Smooth muscle', '20': 'Remove', '21': 'Men Fibro', '22': 'Macro',
        '23': 'Oligo', '24': 'Astro', '25': 'Oligo', '26': 'Remove', '27': 'Oligo',
    }
    
    # Round 3: Refined classification
    round3_annotation = {
        '0': 'Oligo', '1': 'Oligo', '2': 'Astro', '3': 'Oligo', '4': 'OPC',
        '5': 'Micro', '6': 'OPC', '7': 'BEC', '8': 'Astro', '9': 'Fibro',
        '10': 'BEC', '11': 'T cells', '12': 'Oligo', '13': 'Peri', '14': 'OPC',
        '15': 'Smooth Muscle', '16': 'Micro', '17': 'Astro', '18': 'Oligo',
        '19': 'Remove', '20': 'Oligo', '21': 'Oligo', '22': 'OPC',
    }
    
    # Round 4: Final detailed annotation
    final_annotation = {
        '0': 'Oligo', '1': 'Oligo', '2': 'OPC', '3': 'Oligo', '4': 'OPC',
        '5': 'Astro', '6': 'Micro', '7': 'Astro', '8': 'BEC', '9': 'Oligo',
        '10': 'Oligo', '11': 'Astro', '12': 'Fibro', '13': 'BEC', '14': 'Oligo',
        '15': 'OPC', '16': 'T cells', '17': 'Peri', '18': 'Oligo', '19': 'BEC',
        '20': 'SM', '21': 'Micro', '22': 'Micro', '23': 'Oligo', '24': 'Astro',
        '25': 'Astro', '26': 'OPC', '27': 'Remove', '28': 'Oligo', '29': 'Remove',
        '30': 'Oligo', '31': 'Oligo', '32': 'OPC', '33': 'Oligo', '34': 'Oligo',
        '35': 'OPC',
    }
    
    return round1_annotation, round2_annotation, round3_annotation, final_annotation

class NonNeuronalPipeline(SubclusteringPipeline):
    """Specialized pipeline for non-neuronal cells"""
    
    def __init__(self, input_file='adata_nonneuronal_corrected.h5ad'):
        super().__init__(input_file, 'NonNeuronal', 'other')
        
    def run_iterative_analysis(self, metadata_file, max_rounds=4):
        """Run iterative rounds of analysis with progressive refinement"""
        
        print("Starting iterative non-neuronal cell analysis...")
        
        # Get annotation mappings
        round1_map, round2_map, round3_map, final_map = create_nonneuronal_annotation_rounds()
        
        # Round 1: Initial broad classification
        print("\n=== ROUND 1: Initial Broad Classification ===")
        self.load_data()
        self.process_data()
        self.integrate_data()
        self.cluster_data([1.0])
        self.add_metadata(metadata_file)
        
        # Create plots and find markers
        self.create_plots('leiden_res_1_0')
        self.find_marker_genes('leiden_res_1_0')
        
        # Create dotplot
        markers = define_nonneuronal_markers()
        self.create_dotplot(markers, 'leiden_res_1_0')
        
        # Annotate and filter
        self.annotate_clusters(round1_map, 'leiden_res_1_0')
        self.adata.write(f'{self.output_prefix}_round1.h5ad')
        
        # Filter to main cell types
        keep_types = ['Oligo', 'OPC', 'Astro', 'Immune', 'Endo']
        self.filter_cells_by_annotation(keep_types)
        
        # Round 2: Detailed classification
        print("\n=== ROUND 2: Detailed Classification ===")
        self.reprocess_filtered_data()
        self.cluster_data([1.0])
        
        self.create_plots('leiden_res_1_0')
        self.find_marker_genes('leiden_res_1_0')
        self.create_dotplot(markers, 'leiden_res_1_0')
        
        # Annotate with more detail
        self.annotate_clusters(round2_map, 'leiden_res_1_0')
        self.adata.write(f'{self.output_prefix}_round2.h5ad')
        
        # Filter again
        keep_types2 = ['Oligo', 'OPC', 'Astro', 'Micro', 'BEC-cap', 'BEC-arterial', 
                      'Perivasc Fibro', 'T cells', 'Pericyte', 'BEC-venous', 
                      'Smooth muscle', 'Men Fibro', 'Macro']
        self.filter_cells_by_annotation(keep_types2)
        
        # Round 3: Refined classification
        print("\n=== ROUND 3: Refined Classification ===")
        self.reprocess_filtered_data()
        self.cluster_data([1.0])
        
        self.create_plots('leiden_res_1_0')
        self.find_marker_genes('leiden_res_1_0')
        
        # Simplified annotation
        self.annotate_clusters(round3_map, 'leiden_res_1_0')
        self.adata.write(f'{self.output_prefix}_round3.h5ad')
        
        # Filter one more time
        keep_types3 = ['Oligo', 'OPC', 'Astro', 'Micro', 'BEC', 'Fibro', 
                      'T cells', 'Peri', 'Smooth Muscle']
        self.filter_cells_by_annotation(keep_types3)
        
        # Round 4: Final detailed annotation with multiple resolutions
        print("\n=== ROUND 4: Final Multi-Resolution Analysis ===")
        self.reprocess_filtered_data()
        self.cluster_data([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
        
        self.create_plots('leiden_res_1_4')
        self.find_marker_genes('leiden_res_1_2')  # Use intermediate resolution for markers
        
        # Final annotation
        self.annotate_clusters(final_map, 'leiden_res_1_4')
        
        # Create final plots
        self.create_final_plots()
        
        # Save final results  
        self.adata.write(f'{self.output_prefix}_final.h5ad')
        
        print("Non-neuronal cell analysis completed!")
        
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
        self.integrate_data(max_iter=5)  # More iterations for better integration
        
    def create_final_plots(self):
        """Create final publication-quality plots"""
        # Filter to final cell types
        final_types = ['Oligo', 'OPC', 'Astro', 'Micro', 'T cells', 'SM', 'Fibro', 'Peri', 'BEC']
        clean_adata = self.adata[self.adata.obs['cell_class'].isin(final_types)]
        
        # Create clean UMAP with custom colors
        sc.pl.umap(clean_adata, color=['cell_class'], 
                  palette=['#800080', '#9370DB', '#6A5ACD', '#D8BFD8', 
                          '#7E78B8', '#8E89A0', '#BCBDDC', '#8c5c86', '#a64d79'],
                  legend_loc='on data', legend_fontsize=6, frameon=True,
                  save=f'umap_{self.output_prefix}_final_clean.png')
        
        # Create comprehensive dotplot
        marker_genes = {
            'Astro': ['AQP4', 'GFAP'],
            'OPC': ['PDGFRA', 'SOX10', 'BCAN'],
            'Oligo': ['MOG', 'PLP1', 'MOBP', 'OPALIN'],
            'BEC': ['FLT1'],
            'Fibroblast': ['CEMIP'],
            'Pericyte': ['GRM8'],
            'Smooth muscle': ['SLIT3'],
            'Micro': ['APBB1IP'],
            'Tcells': ['PTPRC'],
        }
        
        sc.pl.dotplot(clean_adata, marker_genes, 'cell_class', 
                     dendrogram=True, save=f'dotplot_{self.output_prefix}_final.pdf')
        
        # Create marker expression UMAP
        key_markers = ['AQP4', 'GFAP', 'PDGFRA', 'SOX10', 'MOG', 'PLP1', 
                      'FLT1', 'CEMIP', 'GRM8', 'SLIT3', 'APBB1IP', 'PTPRC']
        sc.pl.umap(clean_adata, color=key_markers, legend_fontsize=6, 
                  frameon=True, save=f'umap_{self.output_prefix}_markers.png')
        
        # Save cluster statistics
        cross_tab = pd.crosstab(clean_adata.obs['cell_class'], clean_adata.obs['orig.ident']).T
        cross_tab.to_csv(f'figures/{self.output_prefix}_final_sample_counts.csv')
        
        cross_tab_region = pd.crosstab(clean_adata.obs['cell_class'], clean_adata.obs['region']).T
        cross_tab_region.to_csv(f'figures/{self.output_prefix}_final_region_counts.csv')
        
    def create_supertype_annotation(self):
        """Create supertype annotations for detailed subclustering"""
        # Map detailed clusters to supertypes for final clean annotation
        supertype_mapping = {
            'Oligo1': 'Oligo', 'Oligo2': 'Oligo', 'Oligo3': 'Oligo', 
            'Oligo4': 'Oligo', 'Oligo5': 'Oligo',
            'OPC1': 'OPC', 'OPC2': 'OPC', 'OPC3': 'OPC',
            'Astro1': 'Astro', 'Astro2': 'Astro', 'Astro3': 'Astro',
            'Micro1': 'Micro', 'Micro2': 'Micro',
            'BEC-venous': 'BEC', 'BEC-arterial': 'BEC', 'BEC-cap': 'BEC',
            'Fibro': 'Fibro', 'Tcell': 'Tcell', 'Peri': 'Peri', 'SM': 'SM',
        }
        return supertype_mapping


def main():
    """Main function for non-neuronal cell analysis"""
    # Initialize pipeline
    pipeline = NonNeuronalPipeline()
    
    # Define metadata file
    metadata_file = '/data/region_samples_list.csv'
    
    # Run complete iterative analysis
    pipeline.run_iterative_analysis(metadata_file)
    
    print("Non-neuronal cell subclustering completed!")


if __name__ == "__main__":
    main()
