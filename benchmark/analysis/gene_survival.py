"""Which panel genes carry survival signal in the bulk cohort, and which cell type expresses them?

For each panel gene present in the bulk cohort:
  1. univariate Cox model  hazard ~ z-scored expression   (lifelines)
  2. Benjamini-Hochberg FDR across all tested genes
  3. assign the gene to the cell type with the highest mean expression in the
     spatial data (data-driven, no hand-curated marker lists); a gene is called
     "immune" if that cell type is an immune population.

coef > 0 -> higher expression = worse survival (risk); coef < 0 -> protective.

Run with the sidish_env conda env. Writes outputs/gene_survival_pdac.csv.
"""
import os
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

ROOT = os.environ.get("SPATIAL_SIDISH_ROOT", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
BULK = os.path.join(ROOT, "PAAD_notebooks/outputs/run_default/processed_bulk.csv")
ADATA = os.path.join(ROOT, "outputs/config2/adata_SIDISH_weibull080.h5ad")
OUT = os.path.join(ROOT, "outputs/gene_survival_pdac.csv")
IMMUNE = {"B Cells", "T Cells", "Macrophages", "CXCL9/10 Cells", "Mast Cells"}

b = pd.read_csv(BULK, index_col=0)
genes = [g for g in b.columns[2:] if b[g].std() > 1e-8]   # drop constant genes
print(f"patients={len(b)}  events={int(b.event.sum())} ({b.event.mean():.1%})  "
      f"genes tested={len(genes)}")

rows = []
for g in genes:
    x = b[g]
    df = pd.DataFrame({"T": b.duration, "E": b.event, "x": (x - x.mean()) / x.std()})
    c = CoxPHFitter().fit(df, "T", "E")
    rows.append((g, c.params_["x"], c.summary.p["x"]))
d = pd.DataFrame(rows, columns=["gene", "coef", "p"]).set_index("gene")
d["q"] = multipletests(d.p, method="fdr_bh")[1]

a = ad.read_h5ad(ADATA)
X = a.X.toarray() if sp.issparse(a.X) else a.X
mean_by_type = pd.DataFrame(X, columns=a.var_names).groupby(
    a.obs["celltype_major"].astype(str).values).mean()
d["top_celltype"] = mean_by_type.idxmax().reindex(d.index)
d["immune"] = d.top_celltype.isin(IMMUNE)
d.sort_values("p").to_csv(OUT)

for label, s in [("nominal p<0.05", d[d.p < 0.05]), ("FDR q<0.05", d[d.q < 0.05])]:
    print(f"\n{label}: {len(s)} genes ({len(s) / len(d):.1%})  "
          f"risk={int((s.coef > 0).sum())} protective={int((s.coef < 0).sum())}  "
          f"immune={int(s.immune.sum())} (protective {int((s.immune & (s.coef < 0)).sum())})")
    print(s.groupby("top_celltype").coef.agg(n="size", risk=lambda c: int((c > 0).sum()))
          .sort_values("n", ascending=False).to_string())
