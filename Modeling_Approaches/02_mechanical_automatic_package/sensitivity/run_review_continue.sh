#!/usr/bin/env bash
# Continue the 2026-10-07 review rerun from a later step (FROM, default analyses): then accept it as the regression reference (the
# configuration changed on purpose, so the strict regression test stops the first pass), rerun the tests,
# run the Hill-kill variant, and rebuild the reports, the findings chapter and the Final Report.
cd "$(dirname "$0")/.." || exit 1
LOG=outputs/pipeline_run_$(date +%Y-%m-%d)_review.log
T=16
say() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }
say "continue from ${FROM:-analyses}"
julia --project=. run_pipeline.jl --from ${FROM:-analyses} --threads $T >> "$LOG" 2>&1
rc=$?; say "pipeline exit $rc"
if [ $rc -ne 0 ]; then
  if grep -q "Stopped: step tests failed" "$LOG"; then
    say "tests failed (expected on a changed configuration: no matching regression reference) — accepting the run"
    julia --project=. test/make_reference.jl >> "$LOG" 2>&1
    julia --project=. run_pipeline.jl --only tests --threads $T >> "$LOG" 2>&1 || { say "tests still failing — stop"; exit 1; }
  else
    say "stopped before the tests — see the log"; exit 1
  fi
fi
say "Hill-kill variant"
bash sensitivity/run_hill_variant.sh >> "$LOG" 2>&1 || say "Hill variant failed (see outputs/sensitivity/hill)"
say "rebuild reports, findings and Final Report"
for step in report14 report_long findings final; do
  julia --project=. run_pipeline.jl --only $step >> "$LOG" 2>&1 || say "$step failed"
done
say "review rerun done"
