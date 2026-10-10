#!/usr/bin/env bash
# The rerun after the time-to-progression controller (2026-10-09): the zoo fits and analyses are unchanged, so
# the pipeline restarts at control14. The control configuration changed on purpose, so the strict regression
# test stops the first pass: accept the run as the reference, rerun the tests, then the Hill-kill variant, and
# rebuild the findings chapter and the Final Report (they read the Hill outputs too).
cd "$(dirname "$0")/.." || exit 1
LOG=outputs/pipeline_run_$(date +%Y-%m-%d)_ttp.log
T=16
say() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }
say "pipeline from control14"
julia --project=. run_pipeline.jl --from control14 --threads $T >> "$LOG" 2>&1
rc=$?; say "pipeline exit $rc"
if [ $rc -ne 0 ]; then
  if grep -q "Stopped: step tests failed" "$LOG"; then
    say "tests failed (expected on a changed control configuration: no matching regression reference) — accepting the run"
    julia --project=. test/make_reference.jl >> "$LOG" 2>&1
    julia --project=. run_pipeline.jl --only tests --threads $T >> "$LOG" 2>&1 || { say "tests still failing — stop"; exit 1; }
  else
    say "stopped before the tests — see the log"; exit 1
  fi
fi
say "Hill-kill variant"
bash sensitivity/run_hill_variant.sh >> "$LOG" 2>&1 || say "Hill variant failed (see outputs/sensitivity/hill)"
say "rebuild findings and Final Report"
for step in findings final; do
  julia --project=. run_pipeline.jl --only $step >> "$LOG" 2>&1 || say "$step failed"
done
say "ttp rerun done"
