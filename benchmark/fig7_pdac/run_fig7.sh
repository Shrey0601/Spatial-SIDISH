#!/bin/bash
# Fig 7 ablations: train the plain VAE (SIDISH phase 1 only) for seeds 0 and 1, then score.
# Run from benchmark/ after the Fig 6 scoring (score7.py reads fig6_pdac/metrics.csv).
#   PY  python with the SIDISH requirements (default: python on PATH)
cd "$(dirname "$0")/.."
PY=${PY:-python}
for s in 0 1; do
  [ -f fig7_pdac/plain_vae_seed$s ] && continue
  PDAC_SEED=$s $PY -u fig7_pdac/train_plain_vae.py 2>&1 | tr '\r' '\n' \
    | grep --line-buffered -vE 'it/s\]|it/s,|^ *$' > fig7_pdac/plain_vae_seed$s.log
done
$PY -u fig7_pdac/score7.py > fig7_pdac/score7.log 2>&1
