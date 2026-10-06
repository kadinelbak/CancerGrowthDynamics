#!/usr/bin/env bash
# Point fits (no bootstrap) + held-out Ce1 prediction for each model-assumption variant, one at a time.
#   bash sensitivity/run_variants.sh [variant ...]      (default: all, in priority order)
# Logs: outputs/sensitivity/<variant>/{fits,heldout}.log; progress lines in outputs/sensitivity/progress.log
cd "$(dirname "$0")/.." || exit 1
VARIANTS=${*:-"V4_ns V1_floor0 V5_ns_sd V6_sd V2_floor25 V3_floor75"}
python sensitivity/make_variants.py > /dev/null || exit 1
P=outputs/sensitivity/progress.log
for v in $VARIANTS; do
  d=outputs/sensitivity/$v
  if [ -f "$d/heldout_ce1.csv" ]; then echo "$(date +%H:%M) $v already done" >> $P; continue; fi
  echo "$(date +%H:%M) $v fits start" >> $P
  STAGES_TOML=$d/stages.toml ZOO_OUT=$d BOOT_B=0 julia -t 20 --project=. run_zoo_fits.jl > $d/fits.log 2>&1
  rc=$?; echo "$(date +%H:%M) $v fits exit $rc" >> $P
  [ $rc -eq 0 ] || continue
  STAGES_TOML=$d/stages.toml ZOO_OUT=$d julia -t 8 --project=. run_heldout.jl > $d/heldout.log 2>&1
  echo "$(date +%H:%M) $v heldout exit $?" >> $P
done
echo "$(date +%H:%M) ALL DONE" >> $P
