"""PDAC arm: headless mirror of config2.ipynb (the run behind 67.8% precision).

Usage:  PDAC_SEED=1 python run_pdac.py pdac_seed1

Exists to answer one question the frozen-trunk sweep cannot: is 67.8% reproducible
end-to-end at a different seed? The sweep found PDAC cell precision spans 4.05-66.87%
across 40 seeds with a frozen trunk, and 67.8% sits at the top of that distribution.
Every LUAD-vs-PDAC conclusion rests on that single number, so it needs replication.

Hyperparameters are copied verbatim from config2.ipynb. SIDISH/ is unmodified.
"""
import sys, os, json, random, time
import numpy as np, pandas as pd, scanpy as sc, torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root, for SIDISH/
from SIDISH import SIDISH as sidish

NAME = sys.argv[1]
SEED = ITE = int(os.environ.get('PDAC_SEED', 0))
OUTPUT_DIR = f'outputs/{NAME}'
FAST = os.environ.get('PDAC_FAST') == '1'

K_NEIGHBORS = 5
LATENT      = 32
LAMBDA_ENV  = 0.1
BETA_EXPR, BETA_NICHE = 1.0, 0.5
EPOCHS_1, I_EPOCHS = (5, 2) if FAST else (225, 20)
EPOCHS_COX  = 50 if FAST else int(os.environ.get('PDAC_COX_EPOCHS', 500))
# S2: the bulk C-index/cell-precision relation is inverse and monotone across 4 seeds,
# so under-fitting the head deliberately is the candidate stabilisation. 500 epochs
# drives train CI to 0.999; 50-100 leaves it near 0.72.
COX_DROPOUT = 0.0
N_ITER, WEIBULL = (1 if FAST else 5), 0.90
TUMOR = ['Tumor Cells', 'CFTR- Tumor Cells']

os.makedirs(OUTPUT_DIR, exist_ok=True)


def set_seed(s):
    np.random.seed(s); torch.manual_seed(s); torch.cuda.manual_seed_all(s); random.seed(s)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


set_seed(SEED)

adata = sc.read_h5ad('outputs/processed_adata_spatial.h5ad')
bulk  = pd.read_csv('outputs/processed_bulk_spatial.csv', index_col=0)
if 'spatial' not in adata.obsm:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values.astype('float32')

chance = 100 * adata.obs.celltype_major.isin(TUMOR).mean()
print(f'=== {NAME} | seed={SEED} ===')
print(f'adata: {adata.n_obs:,} x {adata.n_vars} | bulk {bulk.shape} | chance {chance:.2f}%', flush=True)

t0 = time.time()
sdh = sidish(adata, bulk, 'cpu', seed=ITE, use_spatial_graph=True, k_neighbors=K_NEIGHBORS)
sdh.init_Phase1(epochs=EPOCHS_1, i_epochs=I_EPOCHS, latent_size=LATENT,
                layer_dims=[512, 128], batch_size=512, optimizer='Adam',
                lr=1e-4, lr_3=1e-4, dropout=0, n_covet_genes=50,
                gat_hidden=128, n_heads=4, lambda_env=LAMBDA_ENV,
                beta_expr=BETA_EXPR, beta_niche=BETA_NICHE)
sdh.init_Phase2(EPOCHS_COX, 128, 1e-4, COX_DROPOUT, 0.2, 256)
sdh.train(N_ITER, WEIBULL, 30, OUTPUT_DIR + '/')
runtime = (time.time() - t0) / 60

train_adata = sdh.vae.getEmbedding(clustering=True)
# adata_SIDISH.h5ad is deliberately NOT written: 386 MB/run, regenerable from
# vae_transfer, and disk is tight. Labels + continuous hazard are kept instead --
# the hazard vector is what the cross-seed ensembling needs.
obs = train_adata.obs
try:
    import torch as _t
    from SIDISH.DEEP_COX_ARCHITECTURE import DEEPCOX_ARCHITECTURE
    _cox = DEEPCOX_ARCHITECTURE(128, sdh.vae, COX_DROPOUT)
    _cox.load_state_dict(_t.load(f'{OUTPUT_DIR}/deepCox', map_location='cpu'), strict=False)
    _cox.eval()
    _X = np.asarray(adata.X, dtype=np.float32)
    with _t.no_grad():
        _s = _t.cat([_cox(_t.from_numpy(_X[i:i+8192]).float()).flatten()
                     for i in range(0, _X.shape[0], 8192)]).numpy()
    np.save(f'{OUTPUT_DIR}/cell_hazard.npy', _s)
    print(f'saved cell_hazard.npy  ({_s.shape[0]:,} cells)', flush=True)
except Exception as _e:
    print(f'WARN could not save cell_hazard: {_e}', flush=True)
is_tumor = obs['celltype_major'].isin(TUMOR).astype(int)
is_hr = (obs['SIDISH'] == 'h').astype(int)
tp = int(((is_tumor == 1) & (is_hr == 1)).sum())
fp = int(((is_tumor == 0) & (is_hr == 1)).sum())
fn = int(((is_tumor == 1) & (is_hr == 0)).sum())
precision = 100 * tp / max(tp + fp, 1)
recall = 100 * tp / max(tp + fn, 1)

summary = dict(config=NAME, seed=SEED, cells=int(adata.n_obs), genes=int(adata.n_vars),
               n_high=int(is_hr.sum()), pct_high=round(100 * is_hr.mean(), 2),
               tumor_precision=round(precision, 2), tumor_recall=round(recall, 2),
               chance=round(chance, 2), lift=round(precision / chance, 2),
               test_ci=round(float(sdh.test_ci_Cox[-1]), 4), runtime_min=round(runtime, 1))
json.dump(summary, open(f'{OUTPUT_DIR}/summary.json', 'w'), indent=2)
obs[['celltype_major', 'SIDISH']].to_csv(f'{OUTPUT_DIR}/cell_labels.csv')

print(f'\n===== {NAME} RESULT =====')
for k, v in summary.items():
    print(f'  {k:18s} {v}')
print('\nHigh-risk composition:')
print((obs[obs.SIDISH == 'h'].celltype_major.value_counts(normalize=True) * 100).round(1).to_string())
print('\nconfig2 seed 0 reference: 67.8% precision / 74.5% recall / paper 67.5%', flush=True)
