#!/usr/bin/env bash
# Refit after a change that starts at Stage 3a: the fit step reuses unchanged per-line fits and bootstraps
# (--resume), then every later step recomputes (control checkpoints must not be reused).
cd "$(dirname "$0")/.." || exit 1
LOG=outputs/pipeline_run_$(date +%Y-%m-%d).log
echo "$(date '+%Y-%m-%d %H:%M') refit from s3: fits (resume)" >> "$LOG"
julia --project=. run_pipeline.jl --only fits --resume --threads 20 >> "$LOG" 2>&1 || { echo "$(date '+%H:%M') refit stopped at fits" >> "$LOG"; exit 1; }
echo "$(date '+%Y-%m-%d %H:%M') refit from s3: analyses onward" >> "$LOG"
julia --project=. run_pipeline.jl --from analyses --threads 20 >> "$LOG" 2>&1
echo "$(date '+%Y-%m-%d %H:%M') refit from s3 exit $?" >> "$LOG"
