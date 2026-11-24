#!/bin/bash
# Refit Stage 3a onward with main's Stage 3a candidates (PR #17), recompute control, accept the reference,
# then tests + Final Report. Stops at the first failure.
cd "$(dirname "$0")/.." || exit 1
stop() { echo "[$(date +%H:%M)] CHAIN STOPPED: $1"; exit 1; }
for s in fits analyses; do
  julia -t 12 --project=. run_pipeline.jl --threads 12 --resume --only $s || stop "$s failed"
done
for s in control14 report14 control_long report_long control30 report30 control60 report60 alldoses sweep sweep_all sweep_report findings; do
  REDO=1 CONTROL_FORCE=1 julia -t 12 --project=. run_pipeline.jl --threads 12 --only $s --recompute-control || stop "$s failed"
  if [ $s = control14 ]; then grep -q "checkpoint has Symbol\[\]" outputs/logs/control14.log || stop "control14 loaded an old checkpoint"; echo "control14: empty checkpoint confirmed"; fi
done
echo "[$(date +%H:%M)] make_reference"; julia --project=. test/make_reference.jl || stop "make_reference failed"
julia -t 12 --project=. run_pipeline.jl --threads 12 --from tests || stop "tests/final failed"
echo "[$(date +%H:%M)] CHAIN DONE"
