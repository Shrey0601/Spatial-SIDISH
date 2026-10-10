"""Shared inputs for the Fig-6 benchmark on PDAC Xenium. Replicate r = random 10k-cell
subsample (seed 1000+r), paired with SIDISH full-pipeline seed r. All methods are scored on
the same cells. Scissor/scAB build dense cells x cells matrices, hence the 10k cap."""
import os, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
OUT = 'fig6_pdac/inputs'; os.makedirs(OUT, exist_ok=True)
a = sc.read_h5ad('outputs/processed_adata_spatial.h5ad')
bulk = pd.read_csv('outputs/processed_bulk_spatial.csv', index_col=0)
genes = [g for g in a.var_names if g in bulk.columns]
print(a.shape, 'shared genes', len(genes))
bulk[genes].T.to_csv(f'{OUT}/bulk_expr.csv')                         # genes x patients, log2
bulk[['duration', 'event']].rename(columns={'duration': 'time', 'event': 'status'}).to_csv(f'{OUT}/pheno.csv')
X = a[:, genes].X; X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
for r in range(6):
    idx = np.sort(np.random.default_rng(1000 + r).choice(a.n_obs, 10_000, replace=False))
    np.save(f'{OUT}/idx_{r}.npy', idx)
    pd.DataFrame(X[idx].T.toarray(), index=genes, columns=a.obs_names[idx]).to_csv(f'{OUT}/sc_counts_{r}.csv')
    a.obs.iloc[idx][['x_centroid', 'y_centroid', 'celltype_major']].to_csv(f'{OUT}/sc_meta_{r}.csv')
print('done')
