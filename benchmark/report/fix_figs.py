"""Re-render the three figures that failed the eyeball pass.

fig6  labels overflowed the right edge, x-axis ran to 1e6 (impossible for a P),
      and the P=0.05 tag collided with the frame -> legend instead of inline text,
      axis capped at 1, cohorts ordered with GSE157009 first.
fig5  chance lines were unlabelled.
fig3  the +/-1.3 label sat on top of the chance line.
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, MUTED, SURF = '#0b0b0b', '#52514e', '#8a8984', '#fcfcfb'
FIG = 'report/figs'
plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
    'axes.edgecolor': MUTED, 'axes.linewidth': 0.8, 'axes.labelcolor': INK2,
    'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'font.size': 10,
    'axes.titlesize': 11, 'axes.titleweight': 'bold', 'axes.grid': True,
    'grid.color': '#e6e5e1', 'grid.linewidth': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
})
def save(fig, n):
    fig.tight_layout(); fig.savefig(f'{FIG}/{n}', dpi=200, bbox_inches='tight'); plt.close(fig); print('rewrote', n)

# ---- fig6 ---------------------------------------------------------------
rows = [('GSE157009\n(n=248, squamous)', 5.40e-19, 0.504, 5.24e-05),
        ('GSE37745\n(n=196, mixed NSCLC)', 1.89e-21, 0.489, 0.990)]
fig, ax = plt.subplots(figsize=(7.2, 3.3))
y = np.array([1.0, 0.0])
for i, (lab, claimed, repro, ours) in enumerate(rows):
    yi = y[i]
    ax.scatter(claimed, yi + .20, s=130, color=AQUA, zorder=3, edgecolor=SURF, lw=2)
    ax.scatter(repro,   yi,       s=130, color=ORANGE, marker='s', zorder=3, edgecolor=SURF, lw=2)
    ax.scatter(ours,    yi - .20, s=130, color=BLUE, marker='^', zorder=3, edgecolor=SURF, lw=2)
    ax.annotate(f'{claimed:.0e}', (claimed, yi + .20), xytext=(0, 11), textcoords='offset points',
                ha='center', fontsize=9, color=INK)
    ax.annotate(f'{repro:.2f}', (repro, yi), xytext=(0, 11), textcoords='offset points',
                ha='center', fontsize=9, color=INK)
    ax.annotate(f'{ours:.0e}' if ours < .01 else f'{ours:.2f}', (ours, yi - .20),
                xytext=(0, -20), textcoords='offset points', ha='center', fontsize=9, color=INK)
ax.axvline(0.05, color=MUTED, lw=1.2, ls='--', zorder=1)
ax.text(0.05, 1.46, 'P = 0.05', fontsize=9, color=INK2, ha='center')
ax.set_xscale('log'); ax.set_xlim(1e-24, 3)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.set_ylim(-.55, 1.62)
ax.set_xlabel('log-rank P (lower = stronger stratification)')
ax.set_title("Paper's LUAD external validation does not reproduce")
ax.legend(handles=[
    Line2D([], [], marker='o', ls='', color=AQUA, ms=9, label='paper claims'),
    Line2D([], [], marker='s', ls='', color=ORANGE, ms=9, label="paper's 212 genes, our run"),
    Line2D([], [], marker='^', ls='', color=BLUE, ms=9, label='our SIDISH signature')],
    loc='lower left', bbox_to_anchor=(0.005, -0.02), fontsize=9, labelcolor=INK2, ncol=3)
save(fig, 'fig6_km_reproduction.png')

# ---- fig5 ---------------------------------------------------------------
dp = [0.01, 0.03, 0.12, 0.39, 1.73, 5.58, 17.28, 41.82, 61.85, 67.74]
dl = [4.98, 8.36, 8.82, 11.44, 15.26, 21.34, 24.50, 20.74, 12.95, 6.48]
d = np.arange(10)
fig, ax = plt.subplots(figsize=(6.6, 3.7))
ax.axvspan(8.6, 9.45, color='#f0efe9', zorder=0)
ax.axhline(19.65, color=BLUE, lw=1.1, ls='--', alpha=.5)
ax.axhline(13.54, color=ORANGE, lw=1.1, ls='--', alpha=.5)
ax.text(0.05, 20.5, 'PDAC chance 19.65%', fontsize=8, color=BLUE)
ax.text(0.05, 10.6, 'LUAD chance 13.54%', fontsize=8, color=ORANGE)
ax.plot(d, dp, color=BLUE, lw=2.2, marker='o', ms=6)
ax.plot(d, dl, color=ORANGE, lw=2.2, marker='s', ms=6)
ax.annotate('PDAC — monotone,\npeak at the top decile', (8, dp[8]), xytext=(-12, 6),
            textcoords='offset points', ha='right', color=INK, fontsize=9)
ax.annotate('LUAD — folds,\npeak at decile 6', (6, dl[6]), xytext=(4, 12),
            textcoords='offset points', color=INK, fontsize=9)
ax.text(9.02, 74, 'SIDISH\nreads here', ha='center', fontsize=8.5, color=INK2)
ax.set_xticks(d); ax.set_xlabel('hazard decile (low → high)')
ax.set_ylabel('% tumour cells in decile'); ax.set_ylim(0, 84); ax.set_xlim(-.4, 9.6)
ax.set_title('The hazard axis is monotone in PDAC and folded in LUAD')
save(fig, 'fig5_decile_curves.png')

# ---- fig3 ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), sharey=True)
data = [('PDAC  (chance 19.65%)', [20.7, 12.5, 22.1], [15.1, 2.0, 16.0], BLUE, 19.65,
         'split-driven: 1.12 of variance'),
        ('LUAD  (chance 13.54%)', [15.9, 13.8, 10.2], [16.3, 17.0, 1.3], ORANGE, 13.54,
         'init-driven: 1.09 of variance')]
labs = ['vary both', 'fix split\nvary init', 'vary split\nfix init']
for ax, (title, m, sd, col, ch, note) in zip(axes, data):
    ax.axhline(ch, color=MUTED, lw=1.2, ls='--', zorder=2)
    ax.bar(labs, m, color=col, width=.55, zorder=3)
    ax.errorbar(labs, m, yerr=sd, fmt='none', ecolor=INK2, elinewidth=1.4, capsize=5, zorder=4)
    for i, (mi, s) in enumerate(zip(m, sd)):
        ax.annotate(f'±{s:.1f}', (i, mi + s + 1.6), ha='center', fontsize=9, color=INK)
    ax.text(2.45, ch + 0.8, 'chance', ha='right', fontsize=8.5, color=INK2)
    ax.set_title(title); ax.set_ylim(0, 46)
    ax.text(0.5, -0.32, note, transform=ax.transAxes, ha='center', fontsize=9,
            color=col, fontweight='bold')
axes[0].set_ylabel('tumour-cell precision')
fig.suptitle('Where the seed variance lives — opposite in the two datasets',
             fontsize=11, fontweight='bold', y=1.03)
save(fig, 'fig3_variance_decomposition.png')
