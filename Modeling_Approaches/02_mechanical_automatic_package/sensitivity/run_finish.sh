#!/usr/bin/env bash
# Rebuild the pages after a fit run without refitting: analyses (cached) → both control reports → tests → Final Report.
cd "$(dirname "$0")/.." || exit 1
LOG=outputs/pipeline_run_$(date +%Y-%m-%d).log
for step in analyses report14 report_long tests final; do
  echo "$(date '+%Y-%m-%d %H:%M') --only $step" >> "$LOG"
  julia --project=. run_pipeline.jl --only $step --threads 20 >> "$LOG" 2>&1 || { echo "$(date '+%H:%M') finish stopped at $step" >> "$LOG"; exit 1; }
done
echo "$(date '+%Y-%m-%d %H:%M') finish done" >> "$LOG"
