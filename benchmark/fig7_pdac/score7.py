"""Fig 7 (PDAC Xenium) ablations.
 a-e  SIDISH markers vs bulk-only markers (univariate Cox HR>1, p<0.05 on TCGA; patient
      bootstrap seed r gives 6 replicates). SIDISH rows are read from fig6_pdac/metrics.csv.
 f-j  SIDISH latent (6 seeds) vs simple VAE latent (phase-1 only, no survival): k-means with
      k = #cell types, 20 inits each, on a fixed 20k-cell subsample, scored vs celltype_major.
 k-m  not computed here (needs per-iteration labels/latents).
"""
import os, sys, glob, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, torch, gzip, warnings
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.model_selection import KFold
from sklearn.cluster import KMeans
from sklearn import metrics as skm
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root, for SIDISH/
from SIDISH.VAE import VAE
from SIDISH.Utils import get_spatial_graph_from_adata
PEN = 0.1

bulk = pd.read_csv('outputs/processed_bulk_spatial.csv', index_col=0)
surv = bulk[['duration', 'event']].astype(float); G = bulk.columns[2:]
ext = pd.read_csv('PAAD/GSE224564_expr.txt', sep='\t', index_col=0)
with gzip.open('PAAD/GSE224564_series_matrix.txt.gz', 'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'): titles = [x.strip('"') for x in line.strip().split('\t')[1:]]
        elif '"cancer_death:' in line: ev = [int(x.strip('"').split(': ')[1]) for x in line.strip().split('\t')[1:]]
        elif '"month:' in line: mo = [float(x.strip('"').split(': ')[1]) for x in line.strip().split('\t')[1:]]
es = pd.DataFrame({'duration': mo, 'event': ev}, index=titles); cm = es.index.intersection(ext.columns)
ext_X, es = ext[cm].T, es.loc[cm]


def lr(t, e, hi): return logrank_test(t[hi], t[~hi], e[hi], e[~hi]).p_value


def cox_km(g):   # identical protocol to fig6_pdac/score.py
    d = pd.concat([surv, bulk[g]], axis=1)
    risk = CoxPHFitter(penalizer=PEN).fit(d, 'duration', 'event').predict_partial_hazard(d).values.ravel()
    cv = np.zeros(len(d))
    for tr, te in KFold(5, shuffle=True, random_state=0).split(d):
        cv[te] = CoxPHFitter(penalizer=PEN).fit(d.iloc[tr], 'duration', 'event').predict_partial_hazard(d.iloc[te]).values.ravel()
    ge = [x for x in g if x in ext_X.columns]; s = ext_X[ge].mean(1)
    return dict(n=len(g), km_p=lr(d.duration, d.event, risk > np.median(risk)),
                ci=concordance_index(d.duration, -risk, d.event),
                km_p_cv=lr(d.duration, d.event, cv > np.median(cv)), ci_cv=concordance_index(d.duration, -cv, d.event),
                ext_p=lr(es.duration, es.event, (s >= s.median()).values) if ge else 1.0)


# ---- a-e: bulk-only markers ----
rows = []
for r in range(6):
    bs = surv.sample(frac=1, replace=True, random_state=r)
    bb = bulk.loc[bs.index]
    hits = []
    for g in G:
        try:
            c = CoxPHFitter(penalizer=1e-3).fit(pd.concat([bs.reset_index(drop=True), bb[[g]].reset_index(drop=True)], axis=1),
                                                'duration', 'event')
            if c.hazard_ratios_[g] > 1 and c.summary.p[g] < 0.05: hits.append(g)
        except Exception: pass
    res = cox_km(hits) if hits else dict(n=0, km_p=1, ci=.5, km_p_cv=1, ci_cv=.5, ext_p=1)
    res.update(method='Bulk-only', rep=r); rows.append(res); print(res, flush=True)
f6 = pd.read_csv('fig6_pdac/metrics.csv')
ab = pd.concat([f6[f6.method == 'SIDISH'][list(rows[0].keys())], pd.DataFrame(rows)])
for c in ['km_p', 'km_p_cv', 'ext_p']: ab['neglog10_' + c] = -np.log10(ab[c].clip(1e-300))
ab.to_csv('fig7_pdac/bulk_vs_sidish.csv', index=False)

# ---- f-j: embeddings ----
A = sc.read_h5ad('outputs/processed_adata_spatial.h5ad')
A.X = A.X.toarray() if sp.issparse(A.X) else A.X
A.obsm['spatial'] = A.obs[['x_centroid', 'y_centroid']].values.astype('float32')
vae = VAE(epochs=1, adata=A, z_dim=32, layer_dims=[512, 128], lr=1e-4, dropout=0, device='cpu',
          seed=0, n_covet_genes=50, gat_hidden=128, n_heads=4)
vae.initialize(A, W=np.ones(A.X.shape), batch_size=4096, num_workers=0,
               spatial_graph=get_spatial_graph_from_adata(A, spatial_key='spatial', method='knn', k=5),
               num_neighbors=5, covet_sqrt=np.zeros((A.n_obs, 50, 50), dtype=np.float32))
sub = np.random.default_rng(7).choice(A.n_obs, 20_000, replace=False)
y = A.obs.celltype_major.astype(str).values[sub]; K = len(np.unique(y))
CK = {('SIDISH', 0): 'outputs/config2vae_transfer',
      **{('SIDISH', s): f'outputs/pdac_seed{s}/vae_transfer' for s in range(1, 6)},
      **{('Simple VAE', int(p.split('seed')[-1])): p for p in sorted(p for p in glob.glob('fig7_pdac/plain_vae_seed*') if not p.endswith('.log'))}}
er = []
for (m, s), path in CK.items():
    vae.model.load_state_dict(torch.load(path, map_location='cpu'), strict=False)
    Z = np.asarray(vae.getEmbedding(clustering=True).obsm['latent'])[sub]
    for k in range(20):
        lab = KMeans(K, n_init=1, random_state=k).fit_predict(Z)
        er.append(dict(method=m, seed=s, init=k, ARI=skm.adjusted_rand_score(y, lab),
                       NMI=skm.normalized_mutual_info_score(y, lab), completeness=skm.completeness_score(y, lab),
                       FMI=skm.fowlkes_mallows_score(y, lab),
                       silhouette=skm.silhouette_score(Z, lab, sample_size=5000, random_state=k)))
    print(m, s, pd.DataFrame(er[-20:]).iloc[:, 3:].mean().round(3).to_dict(), flush=True)
E = pd.DataFrame(er); E.to_csv('fig7_pdac/embedding_metrics.csv', index=False)

# ---- figure ----
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
C = {'SIDISH': '#2a78d6', 'Bulk-only': '#8a8a8a', 'Simple VAE': '#eb6834'}
def box(ax, d, title, ms):
    v = [d[m] for m in ms]; b = ax.boxplot(v, patch_artist=True, showfliers=False, medianprops=dict(color='k'))
    for p_, m in zip(b['boxes'], ms): p_.set_facecolor(C[m]); p_.set_alpha(.75)
    p = mannwhitneyu(v[0], v[1], alternative='greater').pvalue if len(v[1]) else np.nan
    ax.set_xticks([1, 2]); ax.set_xticklabels(ms); ax.set_title(f'{title}\nSIDISH > {ms[1]}: P={p:.2g}', fontsize=9)
fig, ax = plt.subplots(2, 5, figsize=(20, 8))
for a_, c, t in zip(ax[0], ['neglog10_km_p', 'ci', 'neglog10_km_p_cv', 'ci_cv', 'neglog10_ext_p'],
                    ['KM -log10P (TCGA in-sample)', 'C-index (in-sample)', 'KM -log10P (TCGA 5-fold CV)',
                     'C-index (5-fold CV)', 'KM -log10P (GSE224564)']):
    box(a_, {m: ab[ab.method == m][c].values for m in ['SIDISH', 'Bulk-only']}, t, ['SIDISH', 'Bulk-only'])
for a_, c in zip(ax[1], ['ARI', 'NMI', 'completeness', 'FMI', 'silhouette']):
    box(a_, {m: E[E.method == m][c].values for m in ['SIDISH', 'Simple VAE']}, c, ['SIDISH', 'Simple VAE'])
fig.suptitle('Fig 7 analogue, PDAC Xenium. Top: SIDISH vs bulk-only markers (6 reps). '
             'Bottom: SIDISH latent (6 seeds) vs simple VAE (phase-1 only), k-means x20 inits', fontsize=11)
plt.tight_layout(rect=(0, 0, 1, .95)); plt.savefig('fig7_pdac/fig7_pdac.png', dpi=140)
print(ab.groupby('method').median(numeric_only=True).round(3).T.to_string())
print(E.groupby('method').median(numeric_only=True).round(3).T.to_string())
