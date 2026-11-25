#!/bin/bash
# Kidney-toxicity side run with extra MPC policies capped at 3 µM (user, 2026-11-22), next to the 1.47 µM MPC:
# 14-day control to outputs/control_kidney_cap3/, 30-day dose + fasting to outputs/control_kidney_cap3_30d/.
cd "$(dirname "$0")/.." || exit 1
stop() { echo "[$(date +%H:%M)] CHAIN STOPPED: $1"; exit 1; }
export CONTROL_TAG=_kidney_cap3 CONTROL_MPC_HIGH_LEVELS=0,0.62,1,1.47,2,3 REDO=1 CONTROL_FORCE=1
mkdir -p outputs/control_kidney_cap3 outputs/control_kidney_cap3_30d
echo "[$(date +%H:%M)] control14 (kidney, cap 3)"; julia -t 12 --project=. run_control.jl > outputs/control_kidney_cap3/run_control.log 2>&1 || stop "run_control failed"
echo "[$(date +%H:%M)] report14"; julia --project=. control_report.jl > outputs/control_kidney_cap3/report.log 2>&1 || stop "control_report failed"
echo "[$(date +%H:%M)] control30 dose + fasting (kidney, cap 3)"; CONTROL_DAYS=30 julia -t 12 --project=. run_control_fasting.jl > outputs/control_kidney_cap3_30d/run_control_fasting.log 2>&1 || stop "run_control_fasting failed"
echo "[$(date +%H:%M)] report30"; CONTROL_DAYS=30 julia --project=. control_fasting_report.jl > outputs/control_kidney_cap3_30d/report.log 2>&1 || stop "control_fasting_report failed"
echo "[$(date +%H:%M)] CHAIN DONE"
