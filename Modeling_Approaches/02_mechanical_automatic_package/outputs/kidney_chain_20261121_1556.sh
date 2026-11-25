#!/bin/bash
# Kidney-toxicity control side run (user, 2026-11-21): the 14-day control (scenario A main, B what-if; 28-day
# look-ahead) to outputs/control_kidney/, then the 30-day dose + fasting run to outputs/control_kidney_30d/.
# outputs/control and the docs/ copies are left as they are (CONTROL_TAG). Stops at the first failure.
cd "$(dirname "$0")/.." || exit 1
stop() { echo "[$(date +%H:%M)] CHAIN STOPPED: $1"; exit 1; }
export CONTROL_TAG=_kidney REDO=1 CONTROL_FORCE=1
echo "[$(date +%H:%M)] control14 (kidney)"; julia -t 12 --project=. run_control.jl > outputs/control_kidney/run_control.log 2>&1 || stop "run_control failed"
echo "[$(date +%H:%M)] report14 (kidney)"; julia --project=. control_report.jl > outputs/control_kidney/report.log 2>&1 || stop "control_report failed"
mkdir -p outputs/control_kidney_30d
echo "[$(date +%H:%M)] control30 dose + fasting (kidney)"; CONTROL_DAYS=30 julia -t 12 --project=. run_control_fasting.jl > outputs/control_kidney_30d/run_control_fasting.log 2>&1 || stop "run_control_fasting failed"
echo "[$(date +%H:%M)] report30 (kidney)"; CONTROL_DAYS=30 julia --project=. control_fasting_report.jl > outputs/control_kidney_30d/report.log 2>&1 || stop "control_fasting_report failed"
echo "[$(date +%H:%M)] CHAIN DONE"
