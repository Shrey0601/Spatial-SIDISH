"""Figures for the SIDISH results report.

Palette: validated categorical slots 1-3 from the dataviz reference instance
(#2a78d6 blue, #eb6834 orange, #1baf7a aqua) -- all-pairs CVD dE 9.2, normal-vision
24.0. Aqua sits below 3:1 on the light surface, so every series is direct-labelled
(the relief rule) rather than identified by colour alone.
"""
import json, os, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#8a8984'
SURF = '#fcfcfb'
FIG = 'report/figs'

plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
    'axes.edgecolor': MUTED, 'axes.linewidth': 0.8, 'axes.labelcolor': INK2,
    'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
    'font.size': 10, 'axes.titlesize': 11, 'axes.titleweight': 'bold',
    'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
})
def save(fig, name):
    fig.tight_layout(); fig.savefig(f'{FIG}/{name}', dpi=200, bbox_inches='tight'); plt.close(fig)
    print('wrote', name)

# ---- 1. PDAC seed spread -------------------------------------------------
seeds, prec, ci = [], [], []
for d, p, c, s in [('config2', 67.80, 0.6676, 0), ('pdac_seed1', 11.74, 0.7716, 1),
                   ('pdac_seed2', 54.39, 0.7224, 2)]:
    seeds.append(s); prec.append(p); ci.append(c)
for f in sorted(glob.glob('outputs/pdac_seed[3-9]/summary.json')):
    j = json.load(open(f)); seeds.append(j['seed']); prec.append(j['tumor_precision']); ci.append(j['test_ci'])

fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.axhline(19.65, color=MUTED, lw=1.4, ls='--', zorder=1)
ax.text(len(seeds) - 0.45, 20.6, 'chance 19.65%', color=INK2, fontsize=9, ha='right')
ax.axhline(67.5, color=AQUA, lw=1.4, ls=':', zorder=1)
ax.text(len(seeds) - 0.45, 68.8, 'paper 67.5%', color=INK2, fontsize=9, ha='right')
x = np.arange(len(seeds))
ax.vlines(x, 0, prec, color='#d6d5d0', lw=1.5, zorder=2)
ax.scatter(x, prec, s=90, color=BLUE, zorder=3, edgecolor=SURF, linewidth=2)
for xi, p in zip(x, prec):
    ax.annotate(f'{p:.1f}%', (xi, p), textcoords='offset points', xytext=(0, 9),
                ha='center', fontsize=9, color=INK)
ax.set_xticks(x); ax.set_xticklabels([f'seed {s}' for s in seeds])
ax.set_ylabel('tumour-cell precision'); ax.set_ylim(0, 78)
ax.set_title('PDAC: identical config, seed varied')
save(fig, 'fig1_pdac_seed_spread.png')

# ---- 2. bulk C-index vs precision ---------------------------------------
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.axhline(19.65, color=MUTED, lw=1.2, ls='--')
ax.text(min(ci) + 0.002, 21, 'chance', color=INK2, fontsize=9)
ax.scatter(ci, prec, s=110, color=BLUE, zorder=3, edgecolor=SURF, linewidth=2)
for c, p, s in zip(ci, prec, seeds):
    ax.annotate(f'seed {s}', (c, p), textcoords='offset points', xytext=(8, -3),
                fontsize=9, color=INK2)
o = np.argsort(ci)
ax.plot(np.array(ci)[o], np.array(prec)[o], color=BLUE, lw=1.6, alpha=.35, zorder=2)
ax.set_xlabel('bulk Cox test C-index'); ax.set_ylabel('tumour-cell precision')
ax.set_title('Better bulk fit, worse cell precision')
save(fig, 'fig2_ci_vs_precision.png')

# ---- 3. variance decomposition ------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), sharey=True)
data = {'PDAC (chance 19.65%)': ([20.7, 12.5, 22.1], [15.1, 2.0, 16.0], BLUE, 19.65),
        'LUAD (chance 13.54%)': ([15.9, 13.8, 10.2], [16.3, 17.0, 1.3], ORANGE, 13.54)}
labs = ['vary both', 'fix split\nvary init', 'vary split\nfix init']
for ax, (title, (m, sd, col, ch)) in zip(axes, data.items()):
    ax.axhline(ch, color=MUTED, lw=1.2, ls='--')
    ax.bar(labs, m, color=col, width=.55, zorder=3)
    ax.errorbar(labs, m, yerr=sd, fmt='none', ecolor=INK2, elinewidth=1.4, capsize=5, zorder=4)
    for i, (mi, s) in enumerate(zip(m, sd)):
        ax.annotate(f'±{s:.1f}', (i, mi + s + 1.2), ha='center', fontsize=9, color=INK)
    ax.set_title(title); ax.set_ylim(0, 45)
axes[0].set_ylabel('tumour-cell precision')
fig.suptitle('Where the seed variance lives — opposite in the two datasets',
             fontsize=11, fontweight='bold', y=1.04)
save(fig, 'fig3_variance_decomposition.png')

# ---- 4. ensemble convergence (two facets, never a dual axis) -------------
k = [1, 2, 4, 5, 8, 10, 20, 40]
pz = [22.3, 21.3, 21.2, 20.8, 20.7, 20.8, 20.5, 20.5]
pr = [22.2, 22.9, 23.1, 23.0, 23.0, 23.3, 23.2, 23.2]
ps = [15.3, 10.0, 7.0, 5.6, 2.9, 4.8, 1.2, np.nan]
lz = [9.8, 10.0, 10.0, 10.1, 10.1, 10.1, 10.1, 10.0]
ls_ = [1.1, 0.6, 0.5, 0.3, 0.4, 0.4, 0.4, np.nan]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 3.4))
a1.axhline(19.65, color=BLUE, lw=1.1, ls='--', alpha=.6)
a1.axhline(13.54, color=ORANGE, lw=1.1, ls='--', alpha=.6)
a1.plot(k, pz, color=BLUE, lw=2, marker='o', ms=6)
a1.plot(k, lz, color=ORANGE, lw=2, marker='s', ms=6)
a1.annotate('PDAC (chance 19.65%)', (k[-1], pz[-1]), xytext=(-4, 9),
            textcoords='offset points', ha='right', color=INK, fontsize=9)
a1.annotate('LUAD (chance 13.54%)', (k[-1], lz[-1]), xytext=(-4, -16),
            textcoords='offset points', ha='right', color=INK, fontsize=9)
a1.set_xscale('log'); a1.set_xticks(k); a1.set_xticklabels(k)
a1.set_xlabel('ensemble size k'); a1.set_ylabel('precision'); a1.set_ylim(0, 28)
a1.set_title('Precision does not improve')
a2.plot(k[:-1], ps[:-1], color=BLUE, lw=2, marker='o', ms=6, label='PDAC')
a2.plot(k[:-1], ls_[:-1], color=ORANGE, lw=2, marker='s', ms=6, label='LUAD')
a2.annotate('PDAC', (k[-2], ps[-2]), xytext=(-2, 9), textcoords='offset points',
            ha='right', color=INK, fontsize=9)
a2.annotate('LUAD', (k[-2], ls_[-2]), xytext=(-2, 7), textcoords='offset points',
            ha='right', color=INK, fontsize=9)
a2.set_xscale('log'); a2.set_xticks(k[:-1]); a2.set_xticklabels(k[:-1])
a2.set_xlabel('ensemble size k'); a2.set_ylabel('sd across disjoint ensembles')
a2.set_title('Variance collapses as 1/√k')
save(fig, 'fig4_ensemble.png')

# ---- 5. decile curves ----------------------------------------------------
dp = [0.01, 0.03, 0.12, 0.39, 1.73, 5.58, 17.28, 41.82, 61.85, 67.74]
dl = [4.98, 8.36, 8.82, 11.44, 15.26, 21.34, 24.50, 20.74, 12.95, 6.48]
fig, ax = plt.subplots(figsize=(6.4, 3.6))
d = np.arange(10)
ax.axhline(19.65, color=BLUE, lw=1.1, ls='--', alpha=.5)
ax.axhline(13.54, color=ORANGE, lw=1.1, ls='--', alpha=.5)
ax.plot(d, dp, color=BLUE, lw=2.2, marker='o', ms=6)
ax.plot(d, dl, color=ORANGE, lw=2.2, marker='s', ms=6)
ax.annotate('PDAC — monotone,\npeak at the top decile', (9, dp[-1]), xytext=(-10, -6),
            textcoords='offset points', ha='right', color=INK, fontsize=9)
ax.annotate('LUAD — folds,\npeak at decile 6', (6, dl[6]), xytext=(6, 10),
            textcoords='offset points', color=INK, fontsize=9)
ax.axvspan(8.6, 9.4, color='#f0efe9', zorder=0)
ax.text(9, 72, 'SIDISH reads here\n(Weibull 0.90)', ha='center', fontsize=8.5, color=INK2)
ax.set_xticks(d); ax.set_xlabel('hazard decile (low → high)')
ax.set_ylabel('tumour cells in decile'); ax.set_ylim(0, 80)
ax.set_title('The hazard axis is monotone in PDAC and folded in LUAD')
save(fig, 'fig5_decile_curves.png')

# ---- 6. KM reproduction --------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 3.2))
rows = [('GSE157009\n(n=248, squamous)', 5.40e-19, 0.504, 5.24e-05),
        ('GSE37745\n(n=196, mixed NSCLC)', 1.89e-21, 0.489, 0.990)]
y = np.arange(len(rows)) * 1.0
for i, (lab, claimed, repro, ours) in enumerate(rows):
    ax.scatter(claimed, y[i] + .18, s=120, color=AQUA, zorder=3, edgecolor=SURF, lw=2)
    ax.scatter(repro, y[i], s=120, color=ORANGE, marker='s', zorder=3, edgecolor=SURF, lw=2)
    ax.scatter(ours, y[i] - .18, s=120, color=BLUE, marker='^', zorder=3, edgecolor=SURF, lw=2)
    ax.annotate(f'{claimed:.1e}  paper claims', (claimed, y[i] + .18), xytext=(9, -3),
                textcoords='offset points', fontsize=8.5, color=INK2)
    ax.annotate(f'{repro:.3f}  their genes, our run', (repro, y[i]), xytext=(9, -3),
                textcoords='offset points', fontsize=8.5, color=INK2)
    ax.annotate(f'{ours:.2e}  our signature', (ours, y[i] - .18), xytext=(9, -3),
                textcoords='offset points', fontsize=8.5, color=INK2)
ax.axvline(0.05, color=MUTED, lw=1.2, ls='--')
ax.text(0.05, 1.62, 'P = 0.05', fontsize=9, color=INK2, ha='center')
ax.set_xscale('log'); ax.set_xlim(1e-23, 1e6)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
ax.set_xlabel('log-rank P'); ax.set_ylim(-.6, 1.8)
ax.set_title("Paper's LUAD external validation does not reproduce")
save(fig, 'fig6_km_reproduction.png')

# ---- 7. LUAD high-risk composition --------------------------------------
comp = [('Endothelial', 31.02, 12.70), ('Alveolar Epithelial', 19.25, 6.07),
        ('Fibroblasts', 12.40, 10.30), ('Macrophages', 10.67, 2.76),
        ('Tumor Cells', 6.35, 13.54), ('Ciliated', 5.28, 1.82),
        ('Monocytes', 3.70, 7.81), ('T Cells', 2.12, 10.89)]
names = [c[0] for c in comp][::-1]
hi = np.array([c[1] for c in comp][::-1]); ov = np.array([c[2] for c in comp][::-1])
yy = np.arange(len(names)); h = .36
fig, ax = plt.subplots(figsize=(6.8, 3.8))
ax.barh(yy + h/2 + .02, hi, height=h, color=ORANGE, zorder=3)
ax.barh(yy - h/2 - .02, ov, height=h, color='#c9c8c2', zorder=3)
for i, (a, b) in enumerate(zip(hi, ov)):
    ax.annotate(f'{a:.1f}%', (a, yy[i] + h/2 + .02), xytext=(4, -3),
                textcoords='offset points', fontsize=8.5, color=INK)
    ax.annotate(f'{b:.1f}%', (b, yy[i] - h/2 - .02), xytext=(4, -3),
                textcoords='offset points', fontsize=8.5, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels(names)
ax.set_xlabel('share of cells'); ax.set_xlim(0, 36)
ax.annotate('high-risk set', (31.0, yy[-1] + h/2 + .02), xytext=(-6, 13),
            textcoords='offset points', ha='right', fontsize=9, color=ORANGE, fontweight='bold')
ax.annotate('whole tissue', (12.7, yy[-1] - h/2 - .02), xytext=(60, -2),
            textcoords='offset points', ha='right', fontsize=9, color=INK2)
ax.set_title('LUAD high-risk cells are vascular and myeloid, not tumour')
save(fig, 'fig7_luad_composition.png')

# ---- 8. LUAD full runs ---------------------------------------------------
fig, ax = plt.subplots(figsize=(6.0, 3.2))
runs = ['seed1\ndp0 ep500', 'config2\ndp0 ep500', 'dp0.0\nep100', 'dp0.5\nep100']
vals = [0.00, 0.65, 2.59, 6.35]
ax.axhline(13.54, color=MUTED, lw=1.4, ls='--')
ax.text(3.45, 14.1, 'chance 13.54%', ha='right', fontsize=9, color=INK2)
ax.axhline(25.02, color=AQUA, lw=1.4, ls=':')
ax.text(3.45, 25.7, 'linear Cox baseline 25.0%', ha='right', fontsize=9, color=INK2)
ax.bar(runs, vals, color=ORANGE, width=.55, zorder=3)
for i, v in enumerate(vals):
    ax.annotate(f'{v:.2f}%', (i, v), xytext=(0, 5), textcoords='offset points',
                ha='center', fontsize=9, color=INK)
ax.set_ylabel('tumour-cell precision'); ax.set_ylim(0, 29)
ax.set_title('LUAD: every full run falls below chance')
save(fig, 'fig8_luad_runs.png')
print('\nall figures written to', FIG)
