#!/bin/bash
# Main control rebuild (user, 2026-11-23): kidney toxicity (Xie et al. 2020) and the MPC at both caps (1.47 and 3 µM).
# Every control step, report and sweep recomputed (fits untouched), then findings, the regression reference, the
# tests and a new Final Report version. Stops at the first failure.
cd "$(dirname "$0")/.." || exit 1
stop() { echo "[$(date +%H:%M)] CHAIN STOPPED: $1"; exit 1; }
for s in control14 report14 control_long report_long control30 report30 control60 report60 sweep sweep_all sweep_report findings; do
  REDO=1 CONTROL_FORCE=1 julia -t 12 --project=. run_pipeline.jl --threads 12 --jobs 1 --only $s --recompute-control || stop "$s failed"
done
echo "[$(date +%H:%M)] make_reference"; julia --project=. test/make_reference.jl || stop "make_reference failed"
julia -t 12 --project=. run_pipeline.jl --threads 12 --jobs 1 --from tests || stop "tests/final failed"
echo "[$(date +%H:%M)] CHAIN DONE"
