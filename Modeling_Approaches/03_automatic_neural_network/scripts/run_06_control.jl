# Stage 5: optimal dose scheduling on the learned models.
# Needs outputs/ude/models.json and outputs/symbolic/symbolic_models.json.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics, JSON3
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"

trajs = N.load_all()
apc = N.area_per_count(trajs)
out = N.outdir("control")
drug = N.make_units(N.training(trajs; role = :drug), apc)

ude_all = N.load_models(joinpath(N.OUT_ROOT, "ude", "models.json"))
sym = JSON3.read(read(joinpath(N.OUT_ROOT, "symbolic", "symbolic_models.json"), String))
optimal_form = sym["optimal"]
N.SYMBOLIC_TERMS[] = (S = Int.(optimal_form.termsS), R = Int.(optimal_form.termsR))
families = [
    "closed form (optimal)" => [N.UDEModel(:symbolic, Float64.(m.theta), Float64.(m.a)) for m in optimal_form.members],
    "UDE (NN kill)" => ude_all["UDE (NN kill)"],
    "Hill (mechanistic)" => ude_all["Hill (mechanistic)"],
]
primary = families[1].second
formula = "kill_S = $(optimal_form.formulaS == "" ? "0" : optimal_form.formulaS); kill_R = $(optimal_form.formulaR == "" ? "0" : optimal_form.formulaR)"

# training support for exposure
vmax_train = maximum(maximum(p[1] for p in N.visited_states(primary[1], drug)[l]; init = 0.0) for l in (-1.0, 1.0))

scenarios = [N.Scenario(name = "A2780 monoculture", S0 = 100.0, R0 = 0.0),
             N.Scenario(name = "50:50 co-culture", S0 = 50.0, R0 = 50.0),
             N.Scenario(name = "25S:75R co-culture", S0 = 25.0, R0 = 75.0)]
settings = [("bounds [0.62, 1.47] (as OptimalControlOne)", 0.62, 1.47),
            ("bounds [0, 1.47] (drug holidays allowed)", 0.0, 1.47)]
n = 13; target = 1.0
lo0, hi0 = 0.62, 1.47
baselines = [
    "constant 1 uM" => fill(1.0, n),
    "front-loaded (max first)" => N.project_mean(collect(range(hi0, lo0; length = n)), target, lo0, hi0),
    "back-loaded (max last)" => N.project_mean(collect(range(lo0, hi0; length = n)), target, lo0, hi0),
    "1 uM on d0-6, off after (as data)" => [i <= 7 ? 1.0 : 0.0 for i in 1:n],
    "off d0-6, 1 uM after (as data)" => [i <= 7 ? 0.0 : 1.0 for i in 1:n],
    "constant 1.47 uM" => fill(1.47, n),
    "no drug" => zeros(n),
]

rows = DataFrame(scenario = String[], setting = String[], schedule = String[], mean_dose = Float64[],
                 total_closed = Float64[], total_closed_sd = Float64[], R_frac_closed = Float64[],
                 total_ude = Float64[], total_hill = Float64[], support = String[])
check = DataFrame(scenario = String[], setting = String[], gradient_obj = Float64[], bbo_obj = Float64[], max_dose_diff = Float64[])
schedules = Dict{Tuple{String, String}, Vector{Float64}}()

function score(sc, doses)
    fin = [N.final_state(m, sc, doses) for m in primary]
    tot = [sum(f) for f in fin]
    res = (mean(tot), std(tot), mean(f[2] / sum(f) for f in fin))
    u = mean(sum(N.final_state(m, sc, doses)) for m in families[2].second)
    h = mean(sum(N.final_state(m, sc, doses)) for m in families[3].second)
    return res..., u, h
end

for sc in scenarios, (sname, lo, hi) in settings
    @info "optimising" sc.name sname
    g = N.optimise_gradient(primary, sc; n, target, lo, hi, starts = FAST ? 3 : 12)
    b = N.optimise_bbo(primary, sc; n, target, lo, hi, max_time = FAST ? 3.0 : 20.0)
    push!(check, (sc.name, sname, g.val, b.val, maximum(abs.(g.sched .- b.sched))))
    best = g.val <= b.val ? g.sched : b.sched
    schedules[(sc.name, sname)] = best
    # schedule optimal under the UDE ensemble, as a cross-model check
    gu = N.optimise_gradient(families[2].second, sc; n, target, lo, hi, starts = FAST ? 2 : 6)
    for (lab, d) in vcat(["OPTIMAL (closed form)" => best, "optimal under UDE" => gu.sched], baselines)
        (lab in first.(baselines) && sname != settings[1][1]) && continue   # baselines once per scenario
        push!(rows, (sc.name, sname, lab, mean(d), score(sc, d)..., N.extrapolation_flags(d; vmax_train)))
    end
end
rows.rel_to_constant = [r.total_closed / rows.total_closed[(rows.scenario .== r.scenario) .& (rows.schedule .== "constant 1 uM")][1] for r in eachrow(rows)]
CSV.write(joinpath(out, "control_results.csv"), rows)
CSV.write(joinpath(out, "optimizer_check.csv"), check)
println(check); println(rows)

# ---------- figures ----------
gr()
ps = []
for sc in scenarios
    p1 = plot(; title = "$(sc.name): schedules", titlefontsize = 8, xlabel = "day", ylabel = "dose (uM)", legend = :outertop, legendfontsize = 6)
    for (sname, lo, hi) in settings
        d = schedules[(sc.name, sname)]
        plot!(p1, 0:n, vcat(d, d[end]); seriestype = :steppost, lw = 2, label = "optimal, $(split(sname, " (")[1])")
    end
    plot!(p1, 0:n, vcat(fill(1.0, n), 1.0); seriestype = :steppost, color = :gray, ls = :dash, label = "constant 1 uM")
    p2 = plot(; title = "$(sc.name): predicted cells (closed form)", titlefontsize = 8, xlabel = "day", ylabel = "cells", yscale = :log10, legend = :outertop, legendfontsize = 6)
    for (lab, d, c) in [("optimal [0.62,1.47]", schedules[(sc.name, settings[1][1])], :royalblue), ("optimal [0,1.47]", schedules[(sc.name, settings[2][1])], :purple), ("constant 1 uM", fill(1.0, n), :gray)]
        for (i, m) in enumerate(primary)
            ts, S, R, _ = N.trajectory(m, sc, d)
            plot!(p2, ts, S .+ R; color = c, alpha = i == 1 ? 1.0 : 0.25, lw = i == 1 ? 2 : 1, label = i == 1 ? "$lab total" : "")
            sc.R0 > 0 && i == 1 && plot!(p2, ts, max.(R, 1e-2); color = c, ls = :dot, label = "$lab resistant")
        end
    end
    push!(ps, p1, p2)
end
savefig(plot(ps...; layout = (length(scenarios), 2), size = (1300, 420 * length(scenarios)), margin = 5Plots.mm), joinpath(out, "control.png"))

N.write_report(joinpath(out, "report.html"), "Stage 5 — optimal dose scheduling",
    [("Set-up", N.para("Same problem as OptimalControlOne: 13 daily doses, mean dose fixed at $(target) uM by the same shift-and-clamp projection, minimise total cells at day 13. Model: recovered closed form ($(formula)), $(length(primary))-member bootstrap ensemble; objective = mean + 0.5·sd of log final total across members (robust). Optimiser: multistart L-BFGS with exact gradients; BlackBoxOptim cross-check. Each schedule is also scored under the UDE and Hill ensembles.")),
     ("Results", N.para("<code>rel_to_constant</code> = final total relative to constant 1 uM (closed form). <code>support</code> flags doses or exposure patterns the training data never covered; treat those predictions as extrapolation.") * N.html_table(rows)),
     ("Optimiser agreement", N.html_table(check)),
     ("Schedules and trajectories", N.figure("control.png", "Optimal schedules and predicted trajectories (thin lines = ensemble members)"))])
println("wrote ", joinpath(out, "report.html"))
