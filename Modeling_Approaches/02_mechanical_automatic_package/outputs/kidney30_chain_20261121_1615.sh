#!/bin/bash
# Kidney-toxicity side run, part 2 (user, 2026-11-21): the 30-day dose + fasting control to outputs/control_kidney_30d/.
cd "$(dirname "$0")/.." || exit 1
stop() { echo "[$(date +%H:%M)] CHAIN STOPPED: $1"; exit 1; }
export CONTROL_TAG=_kidney REDO=1 CONTROL_FORCE=1 CONTROL_DAYS=30
mkdir -p outputs/control_kidney_30d
echo "[$(date +%H:%M)] control30 dose + fasting (kidney)"; julia -t 12 --project=. run_control_fasting.jl > outputs/control_kidney_30d/run_control_fasting.log 2>&1 || stop "run_control_fasting failed"
echo "[$(date +%H:%M)] report30 (kidney)"; julia --project=. control_fasting_report.jl > outputs/control_kidney_30d/report.log 2>&1 || stop "control_fasting_report failed"
echo "[$(date +%H:%M)] CHAIN DONE"
