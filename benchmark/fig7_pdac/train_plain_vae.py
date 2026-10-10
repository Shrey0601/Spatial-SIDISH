"""Fig 7f-j baseline: the 'simple VAE trained without survival guidance'.
This is exactly SIDISH's phase-1 VAE (same config as run_pdac.py / config2), trained on the
PDAC Xenium cells and stopped before any Cox/feedback step. The phase-1 block below mirrors
SIDISH.train() verbatim; SIDISH/ itself is not modified.
Usage: PDAC_SEED=0 python fig7_pdac/train_plain_vae.py
"""
import sys, os, random, time
import numpy as np, pandas as pd, scanpy as sc, torch
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root, for SIDISH/
from SIDISH import SIDISH as sidish
from SIDISH.VAE import VAE
from SIDISH.Utils import compute_covet

SEED = int(os.environ.get('PDAC_SEED', 0))
np.random.seed(SEED); torch.manual_seed(SEED); random.seed(SEED)
adata = sc.read_h5ad('outputs/processed_adata_spatial.h5ad')
bulk = pd.read_csv('outputs/processed_bulk_spatial.csv', index_col=0)
if 'spatial' not in adata.obsm:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values.astype('float32')

t0 = time.time()
s = sidish(adata, bulk, 'cpu', seed=SEED, use_spatial_graph=True, k_neighbors=5)
s.init_Phase1(epochs=int(os.environ.get('PLAIN_EPOCHS', 225)), i_epochs=20, latent_size=32, layer_dims=[512, 128], batch_size=512,
              optimizer='Adam', lr=1e-4, lr_3=1e-4, dropout=0, n_covet_genes=50, gat_hidden=128,
              n_heads=4, lambda_env=0.1, beta_expr=1.0, beta_niche=0.5)
s.init_Phase2(500, 128, 1e-4, 0.0, 0.2, 256)       # as run_pdac.py; keeps RNG state identical up to phase 1
s.num_workers = 8
# ---- mirror of SIDISH.train() phase-1 block (spatial branch) ----
s.covet_sqrt, s.covet_hvg_idx = compute_covet(s.adata, n_hvg=s.n_covet_genes, k=s.k_neighbors,
                                              X_precomputed=getattr(s, 'covet_X_precomputed', None))
s.vae = VAE(epochs=s.epochs_1, adata=s.adata, z_dim=s.latent_size, layer_dims=s.layer_dims, lr=s.lr_1,
            dropout=s.dropout_1, device=s.device, seed=s.seed,
            n_covet_genes=s.n_covet_genes, gat_hidden=s.gat_hidden, n_heads=s.n_heads)
s.vae.initialize(s.adata, W=s.W_matrix, batch_size=s.batch_size, num_workers=s.num_workers,
                 spatial_graph=s.spatial_graph, num_neighbors=s.k_neighbors, covet_sqrt=s.covet_sqrt)
s.vae.train(lambda_env=s.lambda_env, beta_expr=s.beta_expr, beta_niche=s.beta_niche)
os.makedirs('fig7_pdac', exist_ok=True)
OUTP = os.environ.get('PLAIN_OUT', f'fig7_pdac/plain_vae_seed{SEED}'); torch.save(s.vae.model.state_dict(), OUTP)
print(f'saved plain VAE seed {SEED} in {(time.time()-t0)/60:.1f} min', flush=True)
