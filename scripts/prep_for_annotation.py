#!/usr/bin/env python3

import sys
import scanpy as sc
import scanpy.external as sce
import harmonypy as hp
import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp, scipy.io as si
import gzip as gz
import seaborn as sns
from scipy.io import mmwrite

def update_index(row, mappings):
    index = row.name
    new_id = int(row['new_id'])  # Convert new_id from category to int
    if new_id in mappings:
        parts = index.rsplit('-', 1)
        if len(parts) == 2:
            return f'{parts[0]}{mappings[new_id]}'
    return index

# nonneuronal

s = sc.read_h5ad('other_subtype_R4.h5ad') 
m = s.obs
del(s)
c = m.index.values
m['cell_id'] = c
m = pd.DataFrame(m)
m.index = m.index.astype(str)

new_indices = m.apply(lambda row: update_index(row, mappings), axis=1)
m.index = new_indices
c = m.index.values
m['cell_id'] = c

o = sc.read_h5ad('adata_nonneuronal_corrected.h5ad') # subset of original count data 
c = o.obs.index.values
o.obs['cell_id'] = c
o = o[o.obs['cell_id'].isin(m['cell_id'])]
o.obs['cell_class'] = m['cell_class']
o.obs['supertype'] = m['supertype'] 
o.X.shape
o.raw = o

ad.AnnData.write(o, filename = 'other_annotate_corrected.h5ad')

a = o.obs
a = pd.DataFrame(a)
a.to_csv('other_annotate_meta_corrected.csv')

# inhibitory

s = sc.read_h5ad('inn_subtype_R3.h5ad') 
m = s.obs
del(s)
c = m.index.values
m['cell_id'] = c
m = pd.DataFrame(m)
m.index = m.index.astype(str)

new_indices = m.apply(lambda row: update_index(row, mappings), axis=1)
m.index = new_indices
c = m.index.values
m['cell_id'] = c

o = sc.read_h5ad('adata_inn_corrected.h5ad') # subset of original count data 
c = o.obs.index.values
o.obs['cell_id'] = c
o = o[o.obs['cell_id'].isin(m['cell_id'])] 
o.obs['cell_class'] = m['cell_class'] 
o.obs['supertype'] = m['supertype'] 
o.X.shape
o.raw = o

ad.AnnData.write(o, filename = 'inn_annotate_corrected.h5ad')

a = o.obs
a = pd.DataFrame(a)
a.to_csv('inn_annotate_meta_corrected.csv')

# excitatory

s = sc.read_h5ad('exn_subtype_R4.h5ad') 
m = s.obs
del(s)
c = m.index.values
m['cell_id'] = c
m = pd.DataFrame(m)
m.index = m.index.astype(str)

new_indices = m.apply(lambda row: update_index(row, mappings), axis=1)
m.index = new_indices
c = m.index.values
m['cell_id'] = c

o = sc.read_h5ad('adata_exn_corrected.h5ad') # subset of original count data 
c = o.obs.index.values
o.obs['cell_id'] = c
o = o[o.obs['cell_id'].isin(m['cell_id'])]
o.obs['supertype'] = m['supertype'] 
o.obs['cell_class'] = m['cell_class'] 
o.raw = o

u = o[o.obs['cell_class'].isin(['L23IT','L4IT'])]
u.raw.X.shape # (429676, 36591)
l = o[o.obs['cell_class'].isin(['L5IT','L6IT','L5ET','L6CT','L56NP','L6b'])]
l.raw.X.shape # (252328, 36591)

ad.AnnData.write(u, filename = 'exn_upper_annotate_corrected.h5ad')
ad.AnnData.write(l, filename = 'exn_lower_annotate_corrected.h5ad')

um = pd.DataFrame(u.obs)
um.to_csv('exn_upper_annotate_meta_corrected.csv')
lm = pd.DataFrame(l.obs)
lm.to_csv('exn_lower_annotate_meta_corrected.csv')

