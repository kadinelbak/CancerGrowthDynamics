#!/bin/bash
# Control fixes 3 and 4 (user, 2026-11-24), control only on the saved fit set fits_2026-11-10_refit (no fitting):
# 14-day and 30-day dose + fasting control with the progression cost as a step and as a barrier.
# Output folders (CONTROL_TAG; the reports read them): outputs/control_fixes_step, control_fixes_step_30d,
# control_fixes_barrier, control_fixes_barrier_30d. Stops at the first failure.
cd "$(dirname "$0")/.." || exit 1
stop() { echo "[$(date +%H:%M)] CHAIN STOPPED: $1"; exit 1; }
S=$PWD/outputs/fit_sets/fits_2026-11-10_refit
export ZOO_STATE=$S/zoo_fits/zoo_state.jls AN_STATE=$S/zoo_analyses/analyses_state.jls REDO=1 CONTROL_FORCE=1
for cost in step barrier; do
  export CONTROL_TAG=_fixes_$cost CONTROL_PROG_COST=$cost
  mkdir -p outputs/control_fixes_$cost outputs/control_fixes_${cost}_30d
  echo "[$(date +%H:%M)] 14 d $cost"; julia -t 12 --project=. run_control.jl > outputs/control_fixes_$cost/run_control.log 2>&1 || stop "run_control $cost failed"
  julia --project=. control_report.jl > outputs/control_fixes_$cost/report.log 2>&1 || stop "control_report $cost failed"
  echo "[$(date +%H:%M)] 30 d dose + fasting $cost"; CONTROL_DAYS=30 julia -t 12 --project=. run_control_fasting.jl > outputs/control_fixes_${cost}_30d/run_control_fasting.log 2>&1 || stop "run_control_fasting $cost failed"
  CONTROL_DAYS=30 julia --project=. control_fasting_report.jl > outputs/control_fixes_${cost}_30d/report.log 2>&1 || stop "control_fasting_report $cost failed"
done
echo "[$(date +%H:%M)] CHAIN DONE"
