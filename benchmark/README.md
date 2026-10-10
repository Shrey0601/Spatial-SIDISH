# SIDISH on spatial transcriptomics: reproduction and benchmark

Applies **SIDISH** (Jolasun et al., *Nature Communications* 16:11271, 2025) to 10x Xenium data, using the COVET spatial extension in [`../SIDISH`](../SIDISH), and benchmarks it against Scissor, scAB and DEGAS.

Author: Shrey Mehta · Last updated: October 2026

---

## Highlights

### PDAC: the paper's spatial result reproduces

Tumour precision is the share of SIDISH "high-risk" cells that are annotated tumour cells. Chance is 19.65%.

| | Tumour precision | Tumour AUROC |
|---|---|---|
| Paper | 67.5% | — |
| SIDISH, rank-average of 6 full-pipeline runs | **66.8%** | **0.906** |

Averaging the per-cell hazards of several independent runs gives a stable result: every 5-run ensemble lands between 60.1% and 72.0%.

### Cell-level tumour localisation vs other methods (PDAC, 6 replicates)

![Fig 6 benchmark](fig6_pdac/fig6_pdac.png)

| Median | SIDISH | Scissor | scAB | DEGAS |
|---|---|---|---|---|
| Tumour precision of high-risk cells | **56.6%** | 37.3% | 23.4% | 16.9% |
| Per-cell tumour AUROC | **0.87** | 0.54 | 0.57 | 0.48 |
| KM −log10 P, TCGA (paper protocol) | **11.6** | 4.9 | 7.4 | 4.1 |
| C-index, TCGA (paper protocol) | **0.78** | 0.70 | 0.73 | 0.63 |

The KM ranking SIDISH > scAB > Scissor > DEGAS matches the paper's Fig. 6. The cell-level tumour metrics are made possible by Xenium's spatial annotation, and SIDISH leads on them. The full metric table, including cross-validated and external-cohort scores, is in [`fig6_pdac/metrics.csv`](fig6_pdac/metrics.csv). The Fig. 7 ablations are in [`fig7_pdac/`](fig7_pdac).

### High-risk cells mark the tumour niche (PDAC)

High-risk non-tumour cells sit much closer to tumour than low-risk cells of the same type ([`analysis/risk_proximity.py`](analysis/risk_proximity.py)):

| Cell type | Median distance to tumour, high-risk vs low-risk |
|---|---|
| Macrophages | 9.4 vs 23.5 µm |
| Endothelial | 9.9 vs 19.8 µm |
| Fibroblasts | 12.2 vs 29.5 µm |
| T cells | 9.5 vs 46.8 µm |

On 50 µm squares, mean risk rises with tumour fraction (Spearman +0.63).

### Ongoing: LUAD and BRCA

LUAD (Xenium 5K, 275k cells) and BRCA (Xenium HER2+, 159k cells) runs are in progress. In LUAD, high-risk fibroblasts are cancer-associated fibroblasts (POSTN, FAP, COL11A1) sitting at the tumour border: median 25 µm from tumour vs 147 µm for other fibroblasts.

---

## Layout

```
benchmark/
├── run_pdac.py               full SIDISH pipeline, one seed (PDAC_SEED=1 python run_pdac.py pdac_seed1)
├── fig6_pdac/                paper Fig 6: SIDISH vs Scissor / scAB / DEGAS
│   ├── prep.py               6 random 10k-cell subsamples + bulk inputs for the R tools
│   ├── run_tools.R           one tool on one replicate, published defaults
│   ├── run_tools.sh          all tools × replicates, resumable
│   ├── score.py              KM, C-index, external cohort, tumour precision / AUROC
│   ├── plot.py
│   └── metrics.csv, geneset_scores.csv, markers.json, *.png   results
├── fig7_pdac/                paper Fig 7: marker and embedding ablations
│   ├── train_plain_vae.py    SIDISH phase 1 only
│   ├── score7.py
│   ├── run_fig7.sh
│   └── bulk_vs_sidish.csv, embedding_metrics.csv, fig7_pdac.png   results
└── analysis/
    ├── gene_survival.py      per-gene Cox on bulk + which cell type expresses each gene
    └── risk_proximity.py     are high-risk non-tumour cells closer to tumour than their own type?
```

## Running

All scripts run from `benchmark/` and import the package from the repo root. Data is not in the repo; see [`DATA_SETUP.md`](DATA_SETUP.md).

```bash
cd benchmark
pip install -r requirements.txt

# SIDISH full pipeline, one seed (~10 h on CPU)
PDAC_SEED=1 python run_pdac.py pdac_seed1

# Fig 6 benchmark (R env with Scissor, scAB, DEGAS; DEGAS also needs TensorFlow)
python fig6_pdac/prep.py
RSCRIPT=/path/to/Rscript DEGAS_PYTHON=/path/to/python bash fig6_pdac/run_tools.sh
python fig6_pdac/score.py && python fig6_pdac/plot.py

# Fig 7 ablations
bash fig7_pdac/run_fig7.sh

# Neighbourhood check on any run
python analysis/risk_proximity.py --adata <run>/adata_SIDISH.h5ad --celltype-col celltype_major \
    --tumour "Tumor Cells" "CFTR- Tumor Cells" --highrisk-col SIDISH --highrisk-value h
```

## Notes

- The Fig. 6/7 comparison uses the 474-gene Xenium panel and 10k-cell subsamples rather than full-transcriptome scRNA-seq, so absolute values differ from the paper. Compare rankings.
- scAB has 5 of 6 replicates; replicate 1 exceeded the time cap.
