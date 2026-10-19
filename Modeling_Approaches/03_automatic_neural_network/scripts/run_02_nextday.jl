# Stage 1: next-day growth prediction benchmark.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics
const N = NNTrack

t_start = time()
trajs = N.load_all()
res = N.run_nextday(trajs; H = 7)
out = N.outdir("nextday")
CSV.write(joinpath(out, "onestep_metrics.csv"), res.onestep)
CSV.write(joinpath(out, "rollout_metrics.csv"), res.rollout)
CSV.write(joinpath(out, "oof_predictions.csv"), res.oof)
CSV.write(joinpath(out, "onestep_metrics_per_seed.csv"), res.per_seed)
println(res.onestep)

gr()
plts = []
for split in unique(res.rollout.split)
    p = plot(; title = split, xlabel = "days ahead", ylabel = "RMSE of log N", legend = :topleft, titlefontsize = 10)
    for m in unique(res.rollout.model)
        d = res.rollout[(res.rollout.model .== m) .& (res.rollout.split .== split), :]
        plot!(p, d.horizon, d.rmse; yerror = d.rmse_sd_over_seeds, marker = :circle, label = m)
    end
    push!(plts, p)
end
savefig(plot(plts...; layout = (1, length(plts)), size = (1400, 420), margin = 5Plots.mm), joinpath(out, "rollout_rmse.png"))

wide = unstack(res.rollout[res.rollout.horizon .∈ Ref([1, 3, 7]), :], [:model, :split], :horizon, :rmse; renamecols = h -> "h$(h)")
runtime = round((time() - t_start) / 60; digits = 1)
T = "Root-mean-square error of the predicted one-day change in log cell count (0.1 ≈ 10% error per day). Stochastic models (boosted trees, MLP ensemble): mean over $(N.NEXTDAY_SEEDS) seeds, with the SD over seeds."
N.write_report(joinpath(out, "report.html"), "Stage 1 — next-day growth prediction",
    [("What this is", N.para("Target is Δlog N over one day. Features: $(join(string.(N.FEATURES), ", ")). <b>No clock-time features</b>: <code>day</code>, <code>since_switch</code> (time since the last dose change) and <code>cumdose</code> (∫dose dt = dose × day under a constant dose, a treatment clock) were removed on 2026-10-07; the dose history enters only through the filtered exposure v (<span class='eqg-text' data-vars='v,u'>dv/dt = 0.24(u − v)</span>), a saturating state of the input. A row needs measured counts on days t−1, t and t+1 (day 0 is a nominal seeding, not a measurement). Training wells are scored by leave-one-well-out CV; validation wells (ce1 treated co-culture) by models fit on all training wells. Rollouts feed predictions back in for up to 7 days (co-culture partners jointly). Runtime: $(runtime) min.")),
     ("One-step error (log scale)", N.html_table(res.onestep; tip = (T, :lower))),
     ("Multi-step rollout", N.figure("rollout_rmse.png", "Rollout RMSE by horizon (error bars = SD over seeds)"; tip = ("RMSE of log N after feeding predictions back in for 1–7 days.", :lower)) *
        N.html_table(wide; tip = ("Rollout RMSE (log N) at 1, 3 and 7 days ahead, mean over seeds.", :lower)))])
println("wrote ", joinpath(out, "report.html"))
