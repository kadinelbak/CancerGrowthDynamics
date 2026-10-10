#!/usr/bin/env bash
# Finish the time-to-progression rerun after the regression test's NaN fix (2026-10-10): tests on the accepted
# reference, the Hill-kill variant, then the findings chapter and the Final Report.
cd "$(dirname "$0")/.." || exit 1
LOG=outputs/pipeline_run_2026-10-09_ttp.log
T=16
say() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }
say "finish: tests"
julia --project=. run_pipeline.jl --only tests --threads $T >> "$LOG" 2>&1 || { say "tests still failing — stop"; exit 1; }
say "Hill-kill variant"
bash sensitivity/run_hill_variant.sh >> "$LOG" 2>&1 || say "Hill variant failed (see outputs/sensitivity/hill)"
say "rebuild findings and Final Report"
for step in findings final; do
  julia --project=. run_pipeline.jl --only $step >> "$LOG" 2>&1 || say "$step failed"
done
say "ttp rerun done"
