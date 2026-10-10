#!/bin/bash
# Fig 6 benchmark: run Scissor, scAB and DEGAS on replicates 0-5 of the 10k-cell subsample.
# Run from benchmark/ after fig6_pdac/prep.py. Finished outputs are skipped, so it can be resumed.
#   RSCRIPT       Rscript from the R env with Scissor, scAB, DEGAS (default: Rscript on PATH)
#   DEGAS_PYTHON  python with TensorFlow, used by DEGAS (default: python on PATH)
# Time caps: 3 h per run, 1 h for scAB (scAB rep 1 ran past 3 h and was dropped).
cd "$(dirname "$0")/.."
R=${RSCRIPT:-Rscript}
mkdir -p fig6_pdac/results
for r in 0 1 2 3 4 5; do for m in scissor scab degas; do
  [ -f fig6_pdac/results/${m}_${r}.csv ] && continue
  [ $m = scab ] && cap=3600 || cap=10800
  timeout $cap $R fig6_pdac/run_tools.R $m $r > fig6_pdac/results/${m}_${r}.log 2>&1
  grep "^$m $r high" fig6_pdac/results/${m}_${r}.log || echo "$m $r FAILED (cap ${cap}s)"
done; done
