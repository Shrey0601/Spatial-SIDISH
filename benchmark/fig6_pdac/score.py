"""Fig-6 metrics for SIDISH / Scissor / scAB / DEGAS on PDAC Xenium, 6 replicates.

Paper protocol (Methods, 'Kaplan-Meier survival analysis', 'Differential gene expression'):
  markers  = Wilcoxon high-risk vs rest, BH padj < 0.05, log2FC > 1 (paper's standardised cutoff)
  6a  KM   = Cox on markers fit on all TCGA-PAAD patients -> risk score -> median split -> log-rank
  6d  CI   = C-index of that Cox fit (in-sample, as in the paper)
  6g  set  = scanpy score_genes of known cancer genes, evaluated in each method's high-risk cells
Added (not in the paper): 5-fold CV C-index, external KM on GSE224564, cell-level tumour
precision / AUROC against the Xenium annotation.
All methods are scored on the same 10k cells per replicate (fig6_pdac/inputs/idx_r.npy).
"""
import os, json, warnings
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, gzip
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold
warnings.filterwarnings('ignore')

IN, RES = 'fig6_pdac/inputs', 'fig6_pdac/results'
METHODS = ['SIDISH', 'Scissor', 'scAB', 'DEGAS']
SIDISH_DIR = {0: 'outputs/config2', **{s: f'outputs/pdac_seed{s}' for s in range(1, 6)}}
PENALIZER = 0.1          # small ridge so Cox converges with tens of markers on 174 patients

a = sc.read_h5ad('outputs/processed_adata_spatial.h5ad')
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
bulk = pd.read_csv('outputs/processed_bulk_spatial.csv', index_col=0)
surv = bulk[['duration', 'event']].astype(float)
is_tumor = a.obs.celltype_major.isin(['Tumor Cells', 'CFTR- Tumor Cells']).values   # as run_pdac.py

# gene sets for panel g (paper's exact lists are not published)
okb = pd.read_csv(f'{IN}/oncokb_cancerGeneList.tsv', sep='\t')
GS_ONCOKB = sorted(set(okb['Hugo Symbol']) & set(a.var_names))
PDAC_MARKERS = ['KRT19', 'KRT7', 'KRT17', 'KRT8', 'KRT18', 'EPCAM', 'MUC1', 'MUC5AC', 'CEACAM5', 'CEACAM6',
                'MSLN', 'TFF1', 'TFF2', 'S100P', 'LGALS4', 'CLDN18', 'AGR2', 'SOX9', 'ELF3', 'TSPAN8',
                'CDH1', 'MET', 'EGFR', 'ERBB2', 'MYC', 'GATA6', 'KRT5', 'S100A2', 'MKI67', 'TOP2A']
GS_PDAC = [g for g in PDAC_MARKERS if g in a.var_names]
print('gene sets: OncoKB', len(GS_ONCOKB), '| PDAC tumour markers', len(GS_PDAC), GS_PDAC, flush=True)

# external cohort GSE224564 (same loader as gse224564-eval.ipynb)
ext = pd.read_csv('PAAD/GSE224564_expr.txt', sep='\t', index_col=0)
with gzip.open('PAAD/GSE224564_series_matrix.txt.gz', 'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'): titles = [x.strip('"') for x in line.strip().split('\t')[1:]]
        elif '"cancer_death:' in line: ev = [int(x.strip('"').split(': ')[1]) for x in line.strip().split('\t')[1:]]
        elif '"month:' in line: mo = [float(x.strip('"').split(': ')[1]) for x in line.strip().split('\t')[1:]]
ext_surv = pd.DataFrame({'duration': mo, 'event': ev}, index=titles)
common = ext_surv.index.intersection(ext.columns)
ext_X, ext_surv = ext[common].T, ext_surv.loc[common]


def labels(method, r, idx):
    if method == 'SIDISH':
        lab = pd.read_csv(f'{SIDISH_DIR[r]}/cell_labels.csv', index_col=0).SIDISH.values == 'h'
        haz = np.load(f'{SIDISH_DIR[r]}/cell_hazard_rescored.npy')
        return lab[idx], haz[idx]
    f = pd.read_csv(f'{RES}/{method.lower()}_{r}.csv')
    assert (f.cell.values == a.obs_names[idx]).all()
    return f.high.values.astype(bool), f.score.values


def markers(sub, high):
    sub.obs['grp'] = pd.Categorical(np.where(high, 'h', 'b'))
    sc.tl.rank_genes_groups(sub, 'grp', groups=['h'], reference='b', method='wilcoxon')
    rg = sub.uns['rank_genes_groups']
    df = pd.DataFrame({k: rg[k]['h'] for k in ['names', 'logfoldchanges', 'pvals_adj']})
    return df[(df.logfoldchanges > 1) & (df.pvals_adj < 0.05)].names.tolist()


def cox_km(genes):
    g = [x for x in genes if x in bulk.columns]
    if not g: return dict(n=0, km_p=1.0, ci=0.5, ci_cv=0.5, km_p_cv=1.0, ext_p=1.0)
    d = pd.concat([surv, bulk[g]], axis=1)
    cph = CoxPHFitter(penalizer=PENALIZER).fit(d, 'duration', 'event')
    risk = cph.predict_partial_hazard(d).values.ravel()
    hi = risk > np.median(risk)
    km_p = logrank_test(d.duration[hi], d.duration[~hi], d.event[hi], d.event[~hi]).p_value
    ci = concordance_index(d.duration, -risk, d.event)
    cv = np.zeros(len(d))
    for tr, te in KFold(5, shuffle=True, random_state=0).split(d):
        m = CoxPHFitter(penalizer=PENALIZER).fit(d.iloc[tr], 'duration', 'event')
        cv[te] = m.predict_partial_hazard(d.iloc[te]).values.ravel()
    ci_cv = concordance_index(d.duration, -cv, d.event)
    hc = cv > np.median(cv)                              # out-of-fold risk, median split
    km_p_cv = logrank_test(d.duration[hc], d.duration[~hc], d.event[hc], d.event[~hc]).p_value
    ge = [x for x in g if x in ext_X.columns]           # external: mean marker expression, median split
    if ge:
        s = ext_X[ge].mean(1); h2 = s >= s.median()
        ext_p = logrank_test(ext_surv.duration[h2], ext_surv.duration[~h2],
                             ext_surv.event[h2], ext_surv.event[~h2]).p_value
    else: ext_p = 1.0
    return dict(n=len(g), km_p=km_p, ci=ci, ci_cv=ci_cv, km_p_cv=km_p_cv, ext_p=ext_p)


rows, gsrows, mk = [], [], {}
for r in range(6):
    idx = np.load(f'{IN}/idx_{r}.npy')
    sub = a[idx].copy()
    for gs_name, gs in [('OncoKB', GS_ONCOKB), ('PDAC', GS_PDAC)]:
        sc.tl.score_genes(sub, gs, score_name=f'gs_{gs_name}', random_state=0)
    for m in METHODS:
        if m != 'SIDISH' and not os.path.exists(f'{RES}/{m.lower()}_{r}.csv'): continue
        high, score = labels(m, r, idx)
        g = markers(sub, high) if high.sum() >= 10 else []
        mk[f'{m}_{r}'] = g
        res = cox_km(g)
        t = is_tumor[idx]
        res.update(method=m, rep=r, n_high=int(high.sum()),
                   precision=100 * t[high].mean() if high.any() else np.nan,
                   recall=100 * t[high].sum() / t.sum(),
                   auroc=roc_auc_score(t, score) if np.unique(score).size > 1 else np.nan,
                   chance=100 * t.mean())
        rows.append(res)
        for gs_name in ['OncoKB', 'PDAC']:
            v = sub.obs[f'gs_{gs_name}'].values[high]
            gsrows.append(pd.DataFrame({'method': m, 'rep': r, 'set': gs_name, 'score': v}))
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()}, flush=True)

df = pd.DataFrame(rows); df['neglog10_km'] = -np.log10(df.km_p.clip(1e-300))
df['neglog10_km_cv'] = -np.log10(df.km_p_cv.clip(1e-300))
df['neglog10_ext'] = -np.log10(df.ext_p.clip(1e-300))
df.to_csv('fig6_pdac/metrics.csv', index=False)
pd.concat(gsrows).to_csv('fig6_pdac/geneset_scores.csv', index=False)
json.dump(mk, open('fig6_pdac/markers.json', 'w'), indent=1)
pd.set_option('display.width', 200)
print(df.groupby('method')[['n_high', 'n', 'neglog10_km', 'ci', 'neglog10_km_cv', 'ci_cv', 'neglog10_ext', 'precision', 'auroc']]
      .agg(['median', 'mean', 'std']).round(3).T.to_string())
