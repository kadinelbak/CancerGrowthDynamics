# Runs the whole NN-track pipeline in order and logs to outputs/pipeline_run.log.
# Usage (from the 03_automatic_neural_network folder):
#   powershell -NoProfile -ExecutionPolicy Bypass -File scripts/run_all.ps1 [-Fast] [-Out <dir>] [-From <1-9>] [-Threads <n>]
# -From resumes at that step (1 = data ... 8 = learning curve, 9 = Stage 2 hold-out folds, 10 = verdict);
# earlier stages' outputs are reused. The folds run late so Stages 3-7 do not wait for them.
# -Threads caps the Julia threads of every stage (default 6, so a concurrent mechanistic refit keeps its cores).
# Stage 7 (comparison) reads outputs/comparison/mech_heldout_points.csv, written by
#   julia -t 4 --project=../02_mechanical_automatic_package scripts/export_mech_heldout.jl
# which this script runs first when the zoo's zoo_state.jls exists (it only reads the zoo).
param([switch]$Fast, [string]$Out = "", [int]$From = 1, [int]$Threads = 6)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)
if ($Fast) { $env:NN_FAST = "1" } else { Remove-Item Env:NN_FAST -ErrorAction SilentlyContinue }
if ($Out) { $env:NN_OUT = $Out } else { Remove-Item Env:NN_OUT -ErrorAction SilentlyContinue }
$log = if ($Out) { Join-Path $Out "pipeline_run.log" } else { "outputs/pipeline_run.log" }
New-Item -ItemType Directory -Force (Split-Path $log -Parent) | Out-Null
"start $(Get-Date -Format s) from step $From, $Threads threads" | Out-File $log -Append -Encoding ascii
$t = [string]$Threads
$steps = @(
    @("1", "scripts/run_01_data.jl"), @($t, "scripts/run_02_nextday.jl"), @($t, "scripts/run_03_ude.jl"),
    @($t, "scripts/run_04_symbolic.jl"), @($t, "scripts/run_05_importance.jl"), @($t, "scripts/run_06_control.jl"),
    @($t, "scripts/run_07_comparison.jl"), @($t, "scripts/run_08_learning_curve.jl"), @($t, "scripts/run_03b_folds.jl"),
    @("1", "scripts/run_09_verdict.jl"))
for ($i = $From - 1; $i -lt $steps.Count; $i++) {
    $s = $steps[$i]
    if ($i -eq 6 -and -not $Out -and (Test-Path "../02_mechanical_automatic_package/outputs/zoo_fits/zoo_state.jls")) {
        "== export_mech_heldout.jl $(Get-Date -Format s)" | Out-File $log -Append -Encoding ascii
        cmd /c "julia -t 4 --project=../02_mechanical_automatic_package scripts/export_mech_heldout.jl >> `"$log`" 2>&1"
        if ($LASTEXITCODE -ne 0) { "export_mech_heldout.jl failed (exit $LASTEXITCODE); the comparison uses the previous export if any" | Out-File $log -Append -Encoding ascii }
    }
    "== $($s[1]) $(Get-Date -Format s)" | Out-File $log -Append -Encoding ascii
    cmd /c "julia -t $($s[0]) --project=. $($s[1]) >> `"$log`" 2>&1"
    if ($LASTEXITCODE -ne 0) { "FAILED $($s[1]) exit $LASTEXITCODE" | Out-File $log -Append -Encoding ascii; exit 1 }
}
"DONE $(Get-Date -Format s)" | Out-File $log -Append -Encoding ascii
