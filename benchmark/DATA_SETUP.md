# Data setup

No data or model checkpoints are committed. All paths below are relative to `benchmark/`.

## Inputs

| Path | Contents | Used by |
|---|---|---|
| `outputs/processed_adata_spatial.h5ad` | PDAC Xenium, 190,273 cells × 474 genes. Needs `obs['celltype_major']`, `obs['x_centroid']`, `obs['y_centroid']` | everything |
| `outputs/processed_bulk_spatial.csv` | TCGA-PAAD, one row per patient: `duration`, `event`, then the panel genes (log2) | everything |
| `PAAD/GSE224564_expr.txt`, `PAAD/GSE224564_series_matrix.txt.gz` | External PDAC bulk cohort and its survival metadata, from [GEO GSE224564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224564) | `fig6_pdac/score.py`, `fig7_pdac/score7.py` |
| `fig6_pdac/inputs/oncokb_cancerGeneList.tsv` | [OncoKB cancer gene list](https://www.oncokb.org/cancer-genes) | `fig6_pdac/score.py` |

Both processed PDAC files come from `PAAD_notebooks/paad_0_xenium_and_bulk_preprocessing.ipynb`.

## SIDISH run outputs expected by the scoring scripts

`run_pdac.py <name>` writes `outputs/<name>/`. The scoring expects:

| Run | Directory |
|---|---|
| Seed 0 | `outputs/config2/` (the original notebook run) |
| Seeds 1–5 | `outputs/pdac_seed1/` … `outputs/pdac_seed5/` |

| File | Needed by |
|---|---|
| `cell_labels.csv` | `fig6_pdac/score.py` |
| `cell_hazard_rescored.npy` | `fig6_pdac/score.py` |
| `vae_transfer` | `fig7_pdac/score7.py` |
| `summary.json` | `report/make_figs.py` |

Seed 0's checkpoint is read as `outputs/config2vae_transfer`, with no slash, because the notebook saved it under that name.

## Generated (git-ignored)

The following are produced by the scripts and kept out of git:
- `fig6_pdac/inputs/`
- `fig6_pdac/results/`
- `fig7_pdac/plain_vae_seed*`, the plain-VAE weights (about 5 MB each)
- logs
