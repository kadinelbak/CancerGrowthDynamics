#!/usr/bin/env bash
# The whole pipeline (data → fits + bootstrap → analyses → 14-day and long control → tests → Final Report),
# for running detached through run_awake.ps1. Log: outputs/pipeline_run_<date>.log
cd "$(dirname "$0")/.." || exit 1
LOG=outputs/pipeline_run_$(date +%Y-%m-%d).log
echo "$(date '+%Y-%m-%d %H:%M') start" >> "$LOG"
julia --project=. run_pipeline.jl --threads 20 >> "$LOG" 2>&1
echo "$(date '+%Y-%m-%d %H:%M') pipeline exit $?" >> "$LOG"
