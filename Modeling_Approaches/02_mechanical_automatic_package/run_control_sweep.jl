# Control sweep over modelling decisions (PRD "Control barrier and modeling sweep", user 2026-10-27).
# Reruns the normal-medium control study under several burden barriers, thresholds, run lengths, look-ahead
# modes and objectives, so the paper can say which control results come from the biology and which from
# modelling choices. control_sweep_report.jl renders outputs/control_sweep/sweep_outcomes.csv.
#
#   julia -t 12 --project=. run_control_sweep.jl
#   ZOO_STATE=outputs/sensitivity/all_doses/zoo_state.jls SWEEP_DATA="all doses" julia -t 12 --project=. run_control_sweep.jl
#       the same sweep on the fit that keeps the A2780cis 1 µM monoculture (ZOO_NO_EXCLUDE=1 run_zoo_fits.jl);
#       its rows replace only the rows of the same data label, so both data choices end up in one CSV
#
# Variants: each axis changed one at a time from the default (control_config.toml), plus the full
# barrier × run-length cross; every variant at the three nominal seedings and both drug-response scenarios.
# Every run is scored with its own thresholds (B_cancer = B_prog of the variant). Comparators (no treatment,
# constant 1 µM, maximum measured dose, periodic, threshold adaptive) are run under every variant, as the
# threshold policy and every outcome depend on the thresholds and the run length. Data-supported MPC only
# (measured doses, exact search). Days after 14 are beyond the data (column beyond_data).
# Checkpoint: outputs/control_sweep/sweep_state_<data>.jls (REDO=1 recomputes).

include(joinpath(@__DIR__, "control_fasting.jl"))     # includes control_model.jl and control_mpc.jl
using Base.Threads

const SOUT = joinpath(@__DIR__, "outputs", "control_sweep"); mkpath(SOUT)
const DATA_LABEL = get(ENV, "SWEEP_DATA", "without A2780cis 1 µM")
const SKP = joinpath(SOUT, "sweep_state_" * replace(lowercase(DATA_LABEL), r"[^a-z0-9]+" => "_") * ".jls")
const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const LOOK = Int.(CFG["schedule"]["lookahead_blocks"])
const LEVELS = Float64.(CFG["dose"]["levels_data"])
const DATA_DAYS = Float64(CFG["schedule"]["data_days"])
const MPC = "MPC, measured doses only"

# the default variant, from control_config.toml
const D0 = (barrier=String(CFG["cancer"]["barrier"]), threshold=Float64(CFG["cancer"]["B_cancer_fraction_of_K"]),
            days=Float64(CFG["schedule"]["closed_loop_days"]), to_end=Bool(CFG["schedule"]["lookahead_to_end"]),
            objective="time to progression", w=nothing)
const OBJECTIVES = [("time to progression", nothing),
                    ("burden-weighted (qB = qR = 1, no progression term)", [1.0, 1.0, 1.0, 0.1, 0.0, 0.0]),
                    ("no terminal resistant cost (qRT = 0)", [0.1, 0.0, 1.0, 0.1, 10.0, 0.0]),
                    ("3× terminal resistant cost (qRT = 30)", [0.1, 0.0, 1.0, 0.1, 10.0, 30.0])]

function variants()
    V = NamedTuple[]
    add!(axis; kw...) = (v = merge(D0, (axis=axis,), (; kw...)); any(x -> Base.structdiff(x, (axis=0,)) == Base.structdiff(v, (axis=0,)), V) || push!(V, v))
    add!("default")
    for b in CFG["cancer"]["barrier_sweep"], d in CFG["schedule"]["run_days"]   # barrier × run length (full cross)
        add!("barrier × run length"; barrier=String(b), days=Float64(d))
    end
    for t in CFG["cancer"]["B_cancer_sweep"]                                  # threshold (B_cancer and B_prog move together)
        t > 0 && add!("threshold"; threshold=Float64(t))
    end
    for d in CFG["schedule"]["run_days"]                                      # look-ahead: fixed 28 days vs ending with the run
        add!("look-ahead"; days=Float64(d), to_end=!D0.to_end)
    end
    for (o, w) in OBJECTIVES[2:end]
        add!("objective"; objective=o, w=w)
    end
    # the settings before 2026-10-27 (0.8 K, fixed look-ahead); with SWEEP_DATA="all doses" these are the 10 Oct schedules
    add!("pre-2026-10-27 settings"; threshold=0.8, to_end=false)
    V
end
vlabel(v) = "$(v.barrier), $(v.threshold) K, $(round(Int, v.days)) d, $(v.to_end ? "look-ahead to end" : "fixed $(round(Int, sum(LOOK) * Δ0))-day look-ahead"), $(v.objective)"

function run_variant(v, scen, r)
    P = scenario_plant(scen); x0 = x0_for(P, r.N0, r.C0)
    W = default_weights(P; w=v.w, bc_frac=v.threshold, prog_frac=v.threshold, barrier=Symbol(v.barrier))
    n = round(Int, v.days / Δ0)
    pols = Dict{String,Any}()
    for (name, pol) in baseline_policies(P, W)
        pols[name] = run_policy(pol, P, x0, Δ0, n)
    end
    t = time()
    pols[MPC] = closed_loop(P, P, x0, Δ0, n, LOOK, W; mode=:discrete, levels=LEVELS, shrink=v.to_end).U
    secs = time() - t
    (K=P.K, Bc=W.Bc, secs=secs, res=Dict(name => (U=U, o=outcomes(P, x0, U, Δ0, W)) for (name, U) in pols))
end

const V = variants()
const RATIOS = [merge((ratio=r,), nominal_seeding(r)) for r in ("25:75", "50:50", "75:25")]
tasks = [(v, s, r) for v in V for s in (:A, :B) for r in RATIOS]
st = (isfile(SKP) && get(ENV, "REDO", "0") != "1") ? deserialize(SKP) : Dict{Any,Any}()
key(v, s, r) = (vlabel(v), v.axis, s, r.ratio)
todo = [t for t in tasks if !haskey(st, key(t...))]
println("=== control sweep ($DATA_LABEL): $(length(V)) variants × 2 scenarios × 3 ratios = $(length(tasks)) runs, $(length(todo)) to do, $(nthreads()) threads ===")
flush(stdout)
lk = ReentrantLock(); done = Atomic{Int}(0)
@threads :dynamic for i in eachindex(todo)
    v, s, r = todo[i]
    out = run_variant(v, s, r)
    lock(lk) do
        st[key(v, s, r)] = out
        atomic_add!(done, 1)
        @printf("  [%d/%d] %s | %s %s  MPC %.0fs\n", done[], length(todo), vlabel(v), s, r.ratio, out.secs); flush(stdout)
        done[] % 12 == 0 && serialize(SKP, st)
    end
end
serialize(SKP, st)

# ── CSV: the control_outcomes.csv columns plus the variant columns ─────────────
rows = NamedTuple[]
for (v, s, r) in tasks, (pn, x) in sort(collect(st[key(v, s, r)].res); by=first)
    o = x.o; m = st[key(v, s, r)]
    push!(rows, (data=DATA_LABEL, axis=v.axis, variant=vlabel(v), barrier=v.barrier, threshold_fraction_of_K=v.threshold,
                 run_days=v.days, lookahead=v.to_end ? "to end of run" : "fixed", objective=v.objective, beyond_data=v.days > DATA_DAYS,
                 environment="normal", scenario=string(s), ratio=r.ratio, K=m.K, policy=pn, schedule_uM=join(round.(x.U; digits=3), " "),
                 N_final=o.N_final, C_final=o.C_final, total_final=o.total_final, resistant_fraction=o.resistant_fraction,
                 max_total=o.max_total, max_toxicity=o.max_tox, cumulative_dose_uM_days=o.cum_dose,
                 cancer_bound_violated=o.cancer_violation, days_above_B_cancer=o.days_over, toxicity_bound_violated=o.tox_violation,
                 first_progression_day=o.ttp_first, time_to_progression_day=o.ttp, days_progressed=o.days_progressed))
end
new = DataFrame(rows)
f = joinpath(SOUT, "sweep_outcomes.csv")
if isfile(f)                                      # keep the rows of the other data choice
    old = CSV.read(f, DataFrame; stringtype=String)
    old = old[string.(old.data) .!= DATA_LABEL, :]
    new = isempty(old) ? new : vcat(old, new; cols=:union)
end
CSV.write(f, new)
println("Wrote $f ($(nrow(new)) rows)")
