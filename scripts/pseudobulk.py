#!/usr/bin/env python3

import sys
import scanpy as sc
import scanpy.external as sce
import harmonypy as hp
import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp
import scipy.io as si
import gzip as gz
import seaborn as sns
import sc_toolbox as sct
from adpbulk import ADPBulk

# Load new_id data
df = pd.read_csv('/data/region_samples_list.csv', usecols=['sample_id','region','sex','batch','new_id'])
df['sample_id'] = df['sample_id'].str.replace("sample_", "")
df['new_id'] = df['new_id'].str.replace("sample_", "")
id_dict = df.set_index('sample_id')['new_id'].to_dict()

# Load data
other = sc.read_h5ad('other_annotate_corrected.h5ad')
other.obs['orig.ident'] = other.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
print(other.raw.X.shape)
print(other.X[1:100,1:100])
other.obs['new_id'] = other.obs['orig.ident'].map(id_dict).astype('category')
c = pd.crosstab(other.obs['cell_class'], other.obs['new_id']).T
c.to_csv('other_annotate_class_counts_corrected.csv')
c = pd.crosstab(other.obs['supertype'], other.obs['new_id']).T
c.to_csv('other_annotate_supertype_counts_corrected.csv')

inn = sc.read_h5ad('inn_annotate_corrected.h5ad')
inn.obs['orig.ident'] = inn.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
print(inn.X.shape)
print(inn.X[1:100,1:100])
inn.obs['new_id'] = inn.obs['orig.ident'].map(id_dict).astype('category')
c = pd.crosstab(inn.obs['cell_class'], inn.obs['new_id']).T
c.to_csv('inn_annotate_class_counts_corrected.csv')
c = pd.crosstab(inn.obs['supertype'], inn.obs['new_id']).T
c.to_csv('inn_annotate_supertype_counts_corrected.csv')

exn_upper = sc.read_h5ad('exn_upper_annotate_corrected.h5ad')
exn_upper.obs['orig.ident'] = exn_upper.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
print(exn_upper.raw.X.shape)
print(exn_upper.X[1:100,1:100])
exn_upper.obs['new_id'] = exn_upper.obs['orig.ident'].map(id_dict).astype('category')
c = pd.crosstab(exn_upper.obs['cell_class'], exn_upper.obs['new_id']).T
c.to_csv('exn_upper_annotate_class_counts_corrected.csv')
c = pd.crosstab(exn_upper.obs['supertype'], exn_upper.obs['new_id']).T
c.to_csv('exn_upper_annotate_supertype_counts_corrected.csv')

exn_lower = sc.read_h5ad('exn_lower_annotate_corrected.h5ad')
exn_lower.obs['orig.ident'] = exn_lower.obs['orig.ident'].str.replace("_bender_filtered.h5ad", "")
print(exn_lower.raw.X.shape)
print(exn_lower.X[1:100,1:100])
exn_lower.obs['new_id'] = exn_lower.obs['orig.ident'].map(id_dict).astype('category')
c = pd.crosstab(exn_lower.obs['cell_class'], exn_lower.obs['new_id']).T
c.to_csv('exn_lower_annotate_class_counts_corrected.csv')
c = pd.crosstab(exn_lower.obs['supertype'], exn_lower.obs['new_id']).T
c.to_csv('exn_lower_annotate_supertype_counts_corrected.csv')

# Get pseudo-bulk profiles

adpb = ADPBulk(other, ["new_id", "cell_class"], use_raw=True)
other_class_pseudo = adpb.fit_transform()
other_class_pseudo = np.transpose(other_class_pseudo)
other_class_pseudo.to_csv('other_class_pseudo_corrected.csv')

adpb = ADPBulk(other, ["new_id", "supertype"], use_raw=True)
other_class_pseudo = adpb.fit_transform()
other_class_pseudo = np.transpose(other_class_pseudo)
other_class_pseudo.to_csv('other_supertype_pseudo_corrected.csv')

adpb = ADPBulk(inn, ["new_id", "cell_class"], use_raw=True)
inn_class_pseudo = adpb.fit_transform()
inn_class_pseudo = np.transpose(inn_class_pseudo)
inn_class_pseudo.to_csv('inn_class_pseudo_corrected.csv')

adpb = ADPBulk(inn, ["new_id", "supertype"], use_raw=True)
inn_class_pseudo = adpb.fit_transform()
inn_class_pseudo = np.transpose(inn_class_pseudo)
inn_class_pseudo.to_csv('inn_supertype_pseudo_corrected.csv')

adpb = ADPBulk(exn_upper, ["new_id", "cell_class"], use_raw=True)
exn_class_pseudo = adpb.fit_transform()
exn_class_pseudo = np.transpose(exn_class_pseudo)
exn_class_pseudo.to_csv('exn_upper_class_pseudo_corrected.csv')

adpb = ADPBulk(exn_upper, ["new_id", "supertype"], use_raw=True)
exn_class_pseudo = adpb.fit_transform()
exn_class_pseudo = np.transpose(exn_class_pseudo)
exn_class_pseudo.to_csv('exn_upper_supertype_pseudo_corrected.csv')

adpb = ADPBulk(exn_lower, ["new_id", "cell_class"], use_raw=True)
exn_class_pseudo = adpb.fit_transform()
exn_class_pseudo = np.transpose(exn_class_pseudo)
exn_class_pseudo.to_csv('exn_lower_class_pseudo_corrected.csv')

adpb = ADPBulk(exn_lower, ["new_id", "supertype"], use_raw=True)
exn_class_pseudo = adpb.fit_transform()
exn_class_pseudo = np.transpose(exn_class_pseudo)
exn_class_pseudo.to_csv('exn_lower_supertype_pseudo_corrected.csv')
