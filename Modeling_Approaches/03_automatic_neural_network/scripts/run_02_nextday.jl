# Stage 1: next-day growth prediction benchmark.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics
const N = NNTrack

trajs = N.load_all()
res = N.run_nextday(trajs; H = 7)
out = N.outdir("nextday")
CSV.write(joinpath(out, "onestep_metrics.csv"), res.onestep)
CSV.write(joinpath(out, "rollout_metrics.csv"), res.rollout)
CSV.write(joinpath(out, "oof_predictions.csv"), res.oof)
println(res.onestep)

gr()
plts = []
for split in unique(res.rollout.split)
    p = plot(; title = split, xlabel = "days ahead", ylabel = "RMSE of log N", legend = :topleft, titlefontsize = 10)
    for m in unique(res.rollout.model)
        d = res.rollout[(res.rollout.model .== m) .& (res.rollout.split .== split), :]
        plot!(p, d.horizon, d.rmse; marker = :circle, label = m)
    end
    push!(plts, p)
end
savefig(plot(plts...; layout = (1, length(plts)), size = (1400, 420), margin = 5Plots.mm), joinpath(out, "rollout_rmse.png"))

wide = unstack(res.rollout[res.rollout.horizon .∈ Ref([1, 3, 7]), :], [:model, :split], :horizon, :rmse; renamecols = h -> "h$(h)")
N.write_report(joinpath(out, "report.html"), "Stage 1 — next-day growth prediction",
    [("What this is", N.para("Target is Δlog N over one day. Features: $(join(string.(N.FEATURES), ", ")). Training wells are scored by leave-one-well-out CV; validation wells (ce1 treated co-culture) by models fit on all training wells. Rollouts feed predictions back in for up to 7 days (co-culture partners jointly).")),
     ("One-step error (log scale)", N.html_table(res.onestep)),
     ("Multi-step rollout", N.figure("rollout_rmse.png", "Rollout RMSE by horizon") * N.html_table(wide))])
println("wrote ", joinpath(out, "report.html"))
