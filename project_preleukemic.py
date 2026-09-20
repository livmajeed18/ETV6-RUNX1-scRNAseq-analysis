import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import pathlib
import pickle
import utils.proc, utils.plots
import cellproject as cp
from utils.fixes import unclog_umap_caching
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

sc.settings.verbosity = 3
base_figures = './figures/preleukemic_projection/'
base_procdata = './procdata/preleukemic_projection/'
sc.settings.figdir = base_figures
for i in [base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# Load reference
print("Loading reference...")
ref = sc.read('./procdata/04script/combined_filt.h5ad')
with open('./procdata/04script/combined_filt_umapref.pkl', 'rb') as f:
    umapref = pickle.load(f)
ref = ref[ref.obs.data_type == '10x', :].copy()

# Load preleukemic data
print("Loading query...")
query = sc.read('/omics/odcf/analysis/OE0650_projects/rach_seq/src/l172n/obj_mouse_export/obj_mouse.h5ad')

print(f"Reference: {ref.n_obs} cells, {ref.n_vars} genes")
print(f"Query: {query.n_obs} cells, {query.n_vars} genes")

# Unify gene sets
print("Unifying gene sets...")
common_genes = ref.var.index[ref.var.index.isin(query.var_names)]
print(f"Common genes: {len(common_genes)}")
ref_sub = ref[:, common_genes].copy()
ref_sub.X = ref_sub.raw[:, common_genes].X.copy()
del ref_sub.raw
query_sub = query[:, common_genes].copy()

# Combine for batch correction
comb = ref_sub.concatenate(query_sub,
                            batch_key='batch',
                            batch_categories=['ref', 'query'],
                            index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

# Seurat CCA batch correction
comb.obs = comb.obs.astype(str)
print("Running Seurat CCA...")
comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref')

# Project into reference PCA space
print("Projecting...")
sc.pp.scale(ref_sub)
sc.pp.pca(ref_sub, n_comps=50)
sc.pp.neighbors(ref_sub, n_neighbors=15)

query_sub.X = comb_cor[query_sub.obs.index, query_sub.var.index].X.copy()

cp.project_cells(query_sub, ref_sub,
                 obs_columns=['leiden', 'anno_man'],
                 fit_pca=True,
                 scale_data=True)

ref_sub.obsm['X_pca'] = ref_sub.obsm['X_pca_harmony'].copy()

cp.nnregress(query_sub, ref_sub,
             regress=['pca'],
             weighted=True)

cp.project_cells(query_sub, ref_sub,
                 obs_columns=['leiden', 'anno_man'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# Save
query_sub.write(base_procdata + 'preleukemic_projected.h5ad', compression='lzf')
print("Done. preleukemic_projected.h5ad saved.")
