"""Fig-6-style figure for PDAC Xenium. Row 1 = paper's metrics (in-sample TCGA); row 2 = additional
metrics (out-of-fold TCGA, external GSE224564, cell-level ground truth)."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = ['SIDISH', 'Scissor', 'scAB', 'DEGAS']
COL = {'SIDISH': '#2a78d6', 'Scissor': '#eb6834', 'scAB': '#1baf7a', 'DEGAS': '#9b59b6'}
df = pd.read_csv('fig6_pdac/metrics.csv')
gs = pd.read_csv('fig6_pdac/geneset_scores.csv')


def stars(p):
    return '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'n.s.'


def box(ax, data, title, ylab, chance=None):
    vals = [data[m] for m in M]
    b = ax.boxplot(vals, patch_artist=True, widths=.6, showfliers=False, medianprops=dict(color='k'))
    for patch, m in zip(b['boxes'], M): patch.set_facecolor(COL[m]); patch.set_alpha(.75)
    for i, v in enumerate(vals):
        if len(v) <= 50: ax.scatter(np.full(len(v), i + 1) + np.random.default_rng(i).uniform(-.12, .12, len(v)),
                                    v, s=14, c='k', zorder=3)
    top = max(np.nanmax(v) for v in vals); rng_ = top - min(np.nanmin(v) for v in vals) or 1
    for i, m in enumerate(M[1:], start=2):                  # one-sided: SIDISH > other
        p = mannwhitneyu(vals[0], vals[i - 1], alternative='greater').pvalue
        y = top + rng_ * (0.06 + 0.09 * (i - 2))
        ax.plot([1, 1, i, i], [y - rng_*.02, y, y, y - rng_*.02], c='k', lw=.8)
        ax.text((1 + i) / 2, y, stars(p), ha='center', va='bottom', fontsize=8)
    if chance is not None: ax.axhline(chance, ls='--', c='grey', lw=1); ax.text(4.45, chance, 'chance', fontsize=7, va='bottom', ha='right', c='grey')
    ax.set_xticks(range(1, 5)); ax.set_xticklabels(M, fontsize=8); ax.set_title(title, fontsize=10); ax.set_ylabel(ylab, fontsize=9)


by = lambda col: {m: df.loc[df.method == m, col].values for m in M}
fig, ax = plt.subplots(2, 4, figsize=(19, 9))
box(ax[0, 0], by('neglog10_km'), 'a  KM, TCGA in-sample (paper protocol)', '-log10 P')
box(ax[0, 1], by('ci'), 'b  C-index, TCGA in-sample (paper protocol)', 'C-index')
box(ax[0, 2], {m: gs[(gs.method == m) & (gs.set == 'OncoKB')].score.values for m in M},
    'c  OncoKB cancer-gene score in high-risk cells', 'gene set score')
box(ax[0, 3], {m: gs[(gs.method == m) & (gs.set == 'PDAC')].score.values for m in M},
    'd  PDAC tumour-marker score in high-risk cells', 'gene set score')
box(ax[1, 0], by('neglog10_km_cv'), 'e  KM, TCGA out-of-fold (5-fold CV)', '-log10 P')
box(ax[1, 1], by('ci_cv'), 'f  C-index, TCGA out-of-fold (5-fold CV)', 'C-index', chance=.5)
box(ax[1, 2], by('neglog10_ext'), 'g  KM, external GSE224564', '-log10 P', chance=-np.log10(.05))
box(ax[1, 3], by('precision'), 'h  Tumour precision of high-risk cells (Xenium)', '% tumour', chance=df.chance.mean())
fig.text(.005, .74, 'paper\nmetrics', fontsize=11, weight='bold', rotation=90, va='center')
fig.text(.005, .27, 'added\nmetrics', fontsize=11, weight='bold', rotation=90, va='center')
fig.suptitle('PDAC Xenium benchmark, 6 replicates (10k-cell subsample each; SIDISH seed r)  ·  '
             'one-sided Mann-Whitney SIDISH > method', fontsize=11)
plt.tight_layout(rect=(.015, 0, 1, .96)); plt.savefig('fig6_pdac/fig6_pdac.png', dpi=150)

# AUROC as a separate small panel (threshold-free cell-level metric)
fig, a2 = plt.subplots(figsize=(5, 4)); box(a2, by('auroc'), 'Tumour AUROC of per-cell score', 'AUROC', chance=.5)
plt.tight_layout(); plt.savefig('fig6_pdac/fig6_pdac_auroc.png', dpi=150)
print('saved')
