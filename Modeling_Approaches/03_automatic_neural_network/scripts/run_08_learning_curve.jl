# Stage 7: learning curve — does the UDE gain on the mechanistic kill terms as the training data grow?
# Tests the claim "UDEs using NNs were only beneficial when large amounts of data were collected".
#   pool     : the NN track's drug-role training wells (Stage 2: 20k + 30k monocultures, intermittent
#              schedules, low-resource monocultures)
#   subsets  : 25 / 50 / 75 / 100% of the pool, NN_LC_SUBSETS random subsets per fraction (stratified: every
#              subset keeps at least one treated normal-medium well per line); 100% = the whole pool
#   models   : UDE (NN kill; one random seed per subset, NN_SEEDS seeds at 100%), Hill kill, recovered closed
#              form; the competition step always uses all untreated co-culture wells
#   scored on: the same held-out points as Stage 6 — the 30k treated cocultures Ce1 of the zoo's tidy table,
#              days 1-14, from the nominal seeding
# Needs outputs/symbolic/symbolic_models.json. Run: julia -t 6 --project=. scripts/run_08_learning_curve.jl
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics, JSON3, Random, Printf, Serialization
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"
const ITERS = FAST ? (adam_steps = 30, lbfgs_iters = 30) : (adam_steps = 300, lbfgs_iters = 300)
const NSEEDS = FAST ? 2 : parse(Int, get(ENV, "NN_SEEDS", "5"))
const NSUB = FAST ? 2 : parse(Int, get(ENV, "NN_LC_SUBSETS", "4"))   # 4 random subsets per fraction (runtime)
const FRACS = [0.25, 0.5, 0.75, 1.0]
t_start = time()
out = N.outdir("learning_curve")

trajs = N.load_all()
pool = N.make_units(N.training(trajs; role = :drug))
comp = N.make_units(N.training(trajs; role = :competition))
ho = N.make_units(filter(t -> t.source == :ce1, N.load_tidy()))
sym = JSON3.read(read(joinpath(N.OUT_ROOT, "symbolic", "symbolic_models.json"), String))["optimal"]
N.SYMBOLIC_TERMS[] = (S = Int.(sym.termsS), R = Int.(sym.termsR))
coef0 = Float64.(sym.members[1].theta)[N.NB + 1:end]

# stratified random subsets of the pool
line_of(u) = u.members[1].line
okset(us) = all(any(u -> line_of(u) == L && !u.low && any(>(0), u.sched.doses), us) for L in (:S, :R))
rng = N.StableRNG(77)
subsets = NamedTuple[]
for f in FRACS
    k = round(Int, f * length(pool))
    if f == 1.0
        for s in 1:NSEEDS
            push!(subsets, (frac = f, rep = s, seed = s, idx = collect(eachindex(pool))))
        end
    else
        for s in 1:NSUB
            idx = Int[]
            while true
                idx = sort(randperm(rng, length(pool))[1:k])
                okset(pool[idx]) && break
            end
            push!(subsets, (frac = f, rep = s, seed = s, idx))
        end
    end
end
@info "learning curve" pool = length(pool) subsets = length(subsets)

# phase 1: Hill per distinct subset; phase 2: closed form + UDE warm-started from that subset's Hill backbone
# At 100% the data, start values and seeds are exactly Stage 2's main fits (Hill with unit weights = bootstrap
# member 1; UDE seed s = Stage 2 seed s), so those fits are reused from Stage 2's checkpoints when present.
ck(name) = (p = joinpath(N.OUT_ROOT, "ude", "checkpoints", name * ".jls"); isfile(p) ? deserialize(p) : nothing)
hill = Vector{Any}(undef, length(subsets))
Threads.@threads :dynamic for i in eachindex(subsets)
    s = subsets[i]
    # the 100% subsets share the data, and Hill is deterministic: fit it once
    hill[i] = (s.frac == 1.0 && s.rep > 1) ? nothing :
              s.frac == 1.0 && ck("p1_main_hill_1") !== nothing ? ck("p1_main_hill_1") : N.fit_drug_model(:hill, pool[s.idx]; ITERS...)
end
for i in eachindex(subsets)
    hill[i] === nothing && (hill[i] = hill[findfirst(j -> subsets[j].frac == 1.0, eachindex(subsets))])
end
tasks = [(i, kind) for i in eachindex(subsets) for kind in (:symbolic, :nn) if !(kind == :symbolic && subsets[i].frac == 1.0 && subsets[i].rep > 1)]
fits = Dict{Tuple{Int, Symbol}, Any}()
lk = ReentrantLock()
Threads.@threads :dynamic for ti in eachindex(tasks)
    i, kind = tasks[ti]
    s = subsets[i]
    θ0 = kind == :nn ? N.init_theta(:nn; rng = N.StableRNG(1000 + s.seed)) : vcat(zeros(N.NB), coef0)
    θ0[1:N.NB] .= hill[i].θ[1:N.NB]
    reuse = kind == :nn && s.frac == 1.0 ? ck("p2_main_nn_$(s.seed)") : nothing
    f = reuse !== nothing ? reuse : N.fit_drug_model(kind, pool[s.idx]; θinit = θ0, ITERS...)
    m = N.calibrate(kind, f.θ, comp)
    lock(lk) do
        fits[(i, kind)] = m
    end
end
hillm = Vector{Any}(undef, length(subsets))
Threads.@threads :dynamic for i in eachindex(subsets)
    hillm[i] = (subsets[i].frac == 1.0 && subsets[i].rep > 1) ? nothing : N.calibrate(:hill, hill[i].θ, comp)
end
first100 = findfirst(j -> subsets[j].frac == 1.0, eachindex(subsets))
for i in eachindex(subsets)
    hillm[i] === nothing && (hillm[i] = hillm[first100])
    haskey(fits, (i, :symbolic)) || (fits[(i, :symbolic)] = fits[(first100, :symbolic)])
end
@info "fits done" minutes = round((time() - t_start) / 60; digits = 1)

# scores on the held-out points (cells and log), from the nominal seeding
function score(m)
    r = reduce(vcat, [N.forward_residuals(m, u; nominal = true) for u in ho])
    o, p = exp.(r.obs), exp.(r.pred)
    return (rmse = sqrt(mean(abs2, p .- o)), nrmse = sqrt(mean(abs2, p .- o)) / mean(o),
            logrmse = sqrt(mean(abs2, log.(max.(p, 1.0)) .- log.(max.(o, 1.0)))))
end
rows = DataFrame(frac = Float64[], rep = Int[], train_wells = Int[], train_points = Int[], model = String[],
                 rmse_cells = Float64[], nrmse = Float64[], log_rmse = Float64[])
for (i, s) in enumerate(subsets)
    npts = sum(sum(!isnan, y) for u in pool[s.idx] for y in u.y)
    for (name, m) in (("UDE (NN kill)", fits[(i, :nn)]), ("Hill kill (mechanistic)", hillm[i]), ("closed form (mechanistic)", fits[(i, :symbolic)]))
        sc = score(m)
        push!(rows, (s.frac, s.rep, length(s.idx), npts, name, sc.rmse, sc.nrmse, sc.logrmse))
    end
end
CSV.write(joinpath(out, "learning_curve_runs.csv"), rows)

sdz(x) = length(x) > 1 ? std(x) : 0.0
summ = combine(groupby(rows, [:frac, :model]), :train_wells => mean => :train_wells, :rmse_cells => mean => :rmse_cells,
               :rmse_cells => sdz => :rmse_cells_sd, :log_rmse => mean => :log_rmse, :log_rmse => sdz => :log_rmse_sd, nrow => :runs)
# paired differences UDE − mechanistic on the same subset (positive = UDE worse)
pair = DataFrame(frac = Float64[], rep = Int[], vs = String[], d_rmse_cells = Float64[], d_log_rmse = Float64[])
for g in groupby(rows, [:frac, :rep])
    u = g[g.model .== "UDE (NN kill)", :]
    for vs in ("Hill kill (mechanistic)", "closed form (mechanistic)")
        m = g[g.model .== vs, :]
        push!(pair, (g.frac[1], g.rep[1], vs, u.rmse_cells[1] - m.rmse_cells[1], u.log_rmse[1] - m.log_rmse[1]))
    end
end
psum = combine(groupby(pair, [:frac, :vs]), :d_rmse_cells => mean => :d_rmse_cells, :d_rmse_cells => sdz => :d_rmse_cells_sd,
               :d_log_rmse => mean => :d_log_rmse, :d_log_rmse => sdz => :d_log_rmse_sd,
               :d_log_rmse => (x -> count(<(0), x)) => :runs_ude_better_log, :d_rmse_cells => (x -> count(<(0), x)) => :runs_ude_better_cells, nrow => :runs)
CSV.write(joinpath(out, "learning_curve_summary.csv"), summ)
CSV.write(joinpath(out, "learning_curve_paired.csv"), psum)
println(summ); println(psum)

# trend of the paired log difference with training size (least squares slope per 100% of the pool)
trend(vs) = (d = pair[pair.vs .== vs, :]; X = hcat(ones(nrow(d)), d.frac); b = X \ d.d_log_rmse; b[2])
slopes = Dict(vs => trend(vs) for vs in ("Hill kill (mechanistic)", "closed form (mechanistic)"))

gr()
colors = Dict("UDE (NN kill)" => :royalblue, "Hill kill (mechanistic)" => :darkorange, "closed form (mechanistic)" => :crimson)
p1 = plot(; xlabel = "share of training wells used", ylabel = "held-out RMSE (cells)", title = "Held-out Ce1 error vs training size", legend = :topright, titlefontsize = 10)
p2 = plot(; xlabel = "share of training wells used", ylabel = "held-out log-RMSE", title = "Held-out Ce1 log error vs training size", legend = false, titlefontsize = 10)
for (name, c) in colors
    d = sort(summ[summ.model .== name, :], :frac)
    r = rows[rows.model .== name, :]
    jit = name == "UDE (NN kill)" ? -0.01 : name == "Hill kill (mechanistic)" ? 0.0 : 0.01
    scatter!(p1, r.frac .+ jit, r.rmse_cells; color = c, alpha = 0.35, ms = 3, label = "")
    plot!(p1, d.frac, d.rmse_cells; ribbon = d.rmse_cells_sd, color = c, lw = 2, fillalpha = 0.15, marker = :circle, label = name)
    scatter!(p2, r.frac .+ jit, r.log_rmse; color = c, alpha = 0.35, ms = 3, label = "")
    plot!(p2, d.frac, d.log_rmse; ribbon = d.log_rmse_sd, color = c, lw = 2, fillalpha = 0.15, marker = :circle)
end
p3 = plot(; xlabel = "share of training wells used", ylabel = "UDE − mechanistic (log-RMSE)", title = "Paired difference (below 0 = UDE better)", legend = :topright, titlefontsize = 10)
hline!(p3, [0.0]; color = :gray, ls = :dash, label = "")
for vs in ("Hill kill (mechanistic)", "closed form (mechanistic)")
    d = sort(psum[psum.vs .== vs, :], :frac)
    plot!(p3, d.frac, d.d_log_rmse; ribbon = d.d_log_rmse_sd, color = colors[vs], lw = 2, fillalpha = 0.15, marker = :circle, label = "vs $vs")
end
savefig(plot(p1, p2, p3; layout = (1, 3), size = (1500, 440), margin = 5Plots.mm), joinpath(out, "learning_curve.png"))

runtime = round((time() - t_start) / 60; digits = 1)
full_d = psum[psum.frac .== 1.0, :]; low_d = psum[psum.frac .== 0.25, :]
verdict_txt(vs) = begin
    a = low_d[low_d.vs .== vs, :]; b = full_d[full_d.vs .== vs, :]
    "vs $(vs): mean paired log-RMSE difference $(@sprintf("%+.3f", a.d_log_rmse[1])) ± $(@sprintf("%.3f", a.d_log_rmse_sd[1])) at 25% and $(@sprintf("%+.3f", b.d_log_rmse[1])) ± $(@sprintf("%.3f", b.d_log_rmse_sd[1])) at 100% (slope $(@sprintf("%+.3f", slopes[vs])) per 100% of the pool; negative = the UDE gains with more data); UDE better in $(sum(psum[psum.vs .== vs, :runs_ude_better_log])) of $(sum(psum[psum.vs .== vs, :runs])) paired runs (log scale)."
end
N.write_report(joinpath(out, "report.html"), "Learning curve — UDE vs mechanistic kill terms",
    [("Question", N.para("Does the UDE's held-out performance improve <i>relative to</i> the mechanistic kill terms as more training data are used? Training pool: the $(length(pool)) drug-role wells of the NN track (20k and 30k monocultures at 0–1.47 µM, the two intermittent schedules, low-resource monocultures). Each fraction uses $(NSUB) random subsets (stratified so each line keeps a treated normal-medium well); 100% is the whole pool with $(NSEEDS) UDE seeds. The competition step always uses all untreated co-culture wells. Every model is scored on the same $(sum(sum(!isnan, y) for u in ho for y in u.y)) held-out points as Stage 6 (30k Ce1, days 1–14, from the nominal seeding). Runtime: $(runtime) min.")),
     ("Result", N.figure("learning_curve.png", "Held-out error vs share of training wells (points = single runs; line and band = mean ± SD); right: paired UDE − mechanistic difference on the same subsets"; tip = ("Left and middle: held-out error of each model as the training set grows. Right: the UDE's error minus the mechanistic model's error on the same subset; a curve falling below 0 as data grow would support 'the UDE pays off with more data'.", "Lower error is better; in the right panel, below 0 means the UDE is better.")) *
        N.para(join([verdict_txt(vs) for vs in ("Hill kill (mechanistic)", "closed form (mechanistic)")], "<br>")) *
        N.html_table(summ; tip = ("Held-out RMSE (cells) and log-RMSE by training fraction, mean and SD over subsets (and seeds at 100%).", :lower)) *
        N.html_table(psum; tip = ("Paired difference UDE minus mechanistic on the same training subset; runs_ude_better counts the subsets where the UDE had the lower error.", "Lower (more negative) difference favours the UDE; closer to 0 = no difference."))),
     ("Scope", N.para("The largest training set here is still small (≤ $(length(pool)) wells, 4 dose levels, 2 switching schedules), so this curve tests the trend within the data collected, not what would happen with much larger data sets."))])
println("wrote ", joinpath(out, "report.html"), " in $(runtime) min")
