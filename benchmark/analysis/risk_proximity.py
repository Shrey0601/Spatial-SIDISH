"""Do high-risk non-tumour cells sit closer to tumour than low-risk cells of the same type?

The question behind it: when SIDISH flags vessels / macrophages / fibroblasts as
high-risk, are those cells part of the tumour niche, or elsewhere in the tissue?
Comparing high- vs low-risk cells *of the same type* controls for the fact that
some cell types are simply more common near tumour.

Also reports, on square bins, the Spearman correlation between bin mean risk and
bin tumour fraction.

Read-only, no training. Usage (defaults are the PDAC config2 run):
    python analysis/risk_proximity.py
    python analysis/risk_proximity.py --adata <h5ad> --celltype-col <col> \
        --tumour "Tumor" --highrisk-col SIDISH --highrisk-value h \
        --risk-col risk_value --x x_centroid --y y_centroid
"""
import argparse
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr

p = argparse.ArgumentParser()
p.add_argument("--adata", default="outputs/config2/adata_SIDISH_weibull080.h5ad")
p.add_argument("--celltype-col", default="celltype_major")
p.add_argument("--tumour", nargs="+", default=["Tumor Cells", "CFTR- Tumor Cells"])
p.add_argument("--highrisk-col", default="SIDISH")
p.add_argument("--highrisk-value", default="h")
p.add_argument("--risk-col", default="risk_value")
p.add_argument("--x", default="x_centroid")
p.add_argument("--y", default="y_centroid")
p.add_argument("--radius", type=float, default=50.0, help="um, for the 'near tumour' fraction")
p.add_argument("--bin", type=float, default=50.0, help="um, bin size for the bin-level test")
p.add_argument("--min-cells", type=int, default=50, help="skip cell types with fewer high-risk cells")
a = p.parse_args()

obs = ad.read_h5ad(a.adata, backed="r").obs
xy = obs[[a.x, a.y]].to_numpy(float)
ct = obs[a.celltype_col].astype(str).to_numpy()
is_t = np.isin(ct, a.tumour)
is_h = obs[a.highrisk_col].astype(str).to_numpy() == str(a.highrisk_value)
risk = obs[a.risk_col].to_numpy(float)
print(f"cells={len(obs):,}  tumour={is_t.mean():.1%}  high-risk={is_h.mean():.1%}")

# distance from every cell to the nearest tumour cell (tumour cells: nearest *other* tumour cell)
tree = cKDTree(xy[is_t])
d, _ = tree.query(xy, k=2)
dist = np.where(is_t, d[:, 1], d[:, 0])

rows = []
for t in sorted(set(ct) - set(a.tumour)):
    m = ct == t
    h, l = dist[m & is_h], dist[m & ~is_h]
    if len(h) < a.min_cells or len(l) < a.min_cells:
        continue
    rows.append(dict(celltype=t, n_high=len(h), n_low=len(l),
                     med_dist_high=np.median(h), med_dist_low=np.median(l),
                     near_high=(h <= a.radius).mean(), near_low=(l <= a.radius).mean(),
                     p=mannwhitneyu(h, l, alternative="less").pvalue))
res = pd.DataFrame(rows).sort_values("n_high", ascending=False)
res["closer"] = res.med_dist_high < res.med_dist_low
pd.set_option("display.width", 160)
print(f"\nDistance to nearest tumour cell (um), high- vs low-risk cells of the same type"
      f"\n'near' = within {a.radius:g} um; p = one-sided Mann-Whitney (high-risk closer)")
print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# bin level: does bin mean risk track tumour fraction?
b = pd.DataFrame({"bx": (xy[:, 0] // a.bin).astype(int), "by": (xy[:, 1] // a.bin).astype(int),
                  "risk": risk, "t": is_t, "h": is_h})
g = b.groupby(["bx", "by"]).agg(n=("risk", "size"), risk=("risk", "mean"),
                                 tumour=("t", "mean"), frac_h=("h", "mean"))
g = g[g.n >= 5]
rho, pv = spearmanr(g.risk, g.tumour)
print(f"\n{a.bin:g} um bins (>=5 cells): n={len(g):,}  "
      f"Spearman(mean risk, tumour fraction) = {rho:+.3f}  (p={pv:.2g})")
