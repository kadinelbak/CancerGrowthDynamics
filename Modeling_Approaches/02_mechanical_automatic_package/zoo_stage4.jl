# Stage 4: the combination culture — low resource (0.5% FBS) + 1 µM cisplatin, A2780 : A2780cis 1:1,
# 30k (New Datasets/Validation; loaded as `combo` in zoo_common.jl). Fitted as a stage at the user's
# request (2026-09-29); it is not held-out validation.
#
# Every line carries two mechanism modules (zoo_models.jl: modules / line_rhs_multi!):
#   drug module          the Stage 3 carried model of that line (delay chain + Hill kill, crowding, d_N)
#   low-resource module  the Stage 3b carried mechanism of that line (from Stage 1b)
# All inherited values are fixed; candidates add at most a few new parameters:
#   K (shared capacity), χ (drug × low-resource interaction: kill × (1 + χ·u_R/0.955)),
#   ε_N, ε_C (low-resource strength per line), or reset the resistant kill to its Stage 2 value
#   (only offered when the carried Stage 3 model refitted β_C; otherwise it equals the inherited model).
# The model carried forward is chosen by the same rule as every stage (fewest free parameters within
# 2 BIC of the lowest). Include after zoo_common.jl, zoo_models.jl, bootstrap_tools.jl, zoo_selection.jl.

# freeze_slot, freeze_cc: zoo_models.jl. Coculture modifiers (L.cc) of the Stage 3a lines are frozen too.
freeze_line(L, p) = merge(L, (s=NamedTuple{keys(L.s)}(Tuple(freeze_slot(getfield(L.s, k), p) for k in keys(L.s))),),
                          cc_on(L) ? (cc=freeze_cc(L.cc, p),) : NamedTuple())

# (name, label, crowding source, resistant kill on, free parameters) from the combination stage of stages.toml
const S4_SPECS = let st = only(s for s in TOML.parsefile(STAGES_FILE)["stage"] if s["kind"] == "combination")
    [(c["name"], c["label"], c["crowd"] == "s3b" ? :S3b : :S3, get(c, "reskill", false), Symbol[Symbol(x) for x in c["free"]]) for c in st["candidate"]]
end
# stage-level Norton–Simon coupling of the drug kill to growth (stages.toml kill_growth_coupled)
const S4_NS = get(only(s for s in TOML.parsefile(STAGES_FILE)["stage"] if s["kind"] == "combination"), "kill_growth_coupled", false)
const S4_PNAME = Dict(:K => "K", :χ => "χ", :χN => "χ_N", :χC => "χ_C", :εN => "ε_N^R", :εC => "ε_C^R")

combo_datasets(obs) = [dataset(keep(series(COMBO_SOURCES.N)), obs[1]), dataset(keep(series(COMBO_SOURCES.C)), obs[2])]

# The data start at day 0 (the seeding, 50 cells per line); exposure starts then, so the hidden drug and
# stress states start at 0 and the starting counts follow the seeding rule (zoo_data.jl x0_info): nominal,
# or estimated where day 1 is more than 3-fold off it. (This replaces the earlier start at day 1 with
# hidden states simulated from the candidates' starting guesses.)
combo_x0() = [x0_info(series(COMBO_SOURCES.N)), x0_info(series(COMBO_SOURCES.C))]

function stage4_cands(c3, c3b, βC2)
    LN3, LC3 = freeze_line(c3.sys.lines[1], c3.params), freeze_line(c3.sys.lines[2], c3.params)
    MN, MC = freeze_line(c3b.sys.lines[1], c3b.params), freeze_line(c3b.sys.lines[2], c3b.params)
    a3 = [[freeze_slot(a, c3.params) for a in row] for row in c3.sys.alpha]
    a3b = [[freeze_slot(a, c3b.params) for a in row] for row in c3b.sys.alpha]
    K3, K3b = LN3.s.K.v, freeze_slot(c3b.sys.lines[1].s.K, c3b.params).v
    X0 = [series(COMBO_SOURCES.N).nom, series(COMBO_SOURCES.C).nom]; X0info = combo_x0()
    ymax = maximum(vcat(combo.yN, combo.yC))
    cands = NamedTuple[]
    βC3 = LC3.s.β.v
    for (name, label, crowd, reskill, free) in S4_SPECS
        reskill && βC3 == βC2 && continue           # identical to the inherited candidate
        idx = Dict(k => i for (i, k) in enumerate(free))
        sl(k, v) = haskey(idx, k) ? Slot(idx[k], 0.0) : Slot(0, v)
        Kv = crowd == :S3 ? K3 : K3b
        LN = merge(LN3, (s=merge(LN3.s, (K=sl(:K, Kv),)),))
        LC = merge(LC3, (s=merge(LC3.s, reskill ? (K=sl(:K, Kv), β=Slot(0, βC2)) : (K=sl(:K, Kv),)),))
        mN = merge(MN, (s=merge(MN.s, (ε=sl(:εN, MN.s.ε.v),)),))
        mC = merge(MC, (s=merge(MC.s, (ε=sl(:εC, MC.s.ε.v),)),))
        χN = haskey(idx, :χ) ? Slot(idx[:χ], 0.0) : sl(:χN, 0.0)
        χC = haskey(idx, :χ) ? Slot(idx[:χ], 0.0) : sl(:χC, 0.0)
        lines = [merge(LN, (extra=Any[mN], chi=χN)), merge(LC, (extra=Any[mC], chi=χC))]
        S4_NS && (lines = [merge(l, (ns=true,)) for l in lines])
        alpha = crowd == :S3 ? a3 : a3b
        sys = build_system(lines, [(D=COMBO_DOSE, R=U_R, lines=[1, 2])], alpha)
        model = make_system_model(sys)
        bounds = [k == :K ? (0.5ymax, 20ymax) : k in (:χ, :χN, :χC) ? (-0.95, 10.0) : (0.0, 20.0) for k in free]
        start(k) = k == :K ? Kv : k in (:εN, :εC) ? 1.0 : 0.0
        p0s = [[start(k) for k in free], [k == :K ? 1.5Kv : k in (:εN, :εC) ? 0.5 : 1.0 for k in free],
               [k == :K ? 0.7Kv : k in (:εN, :εC) ? 2.0 : -0.5 for k in free]]
        p0s = [[clamp(x, lo + 1e-3 * (hi - lo), hi - 1e-3 * (hi - lo)) for (x, (lo, hi)) in zip(p, bounds)] for p in p0s]
        u0 = sys_u0(sys, [X0])
        push!(cands, (name=name, label=label, family=isempty(free) ? "Combination, inherited" : "Combination, refitted", sys=sys, model=model,
                      u0=u0, pnames=[S4_PNAME[k] for k in free], bounds=bounds, p0s=p0s, slots=free, crowd=crowd, reskill=reskill, βC2=βC2,
                      x0_fields(sys, X0info)...))
    end
    cands
end

function fit_stage4(cands, ds)
    rows = Vector{Any}(undef, length(cands))
    Threads.@threads :dynamic for i in eachindex(cands)
        rows[i] = fit_cand(cands[i], ds)                # bootstrap_tools.jl (estimated starting values included)
    end
    order_rule(sort!([r for r in rows if isfinite(r.bic)]; by=r -> r.bic))
end

function stage4_group(c3, c3b, βC2)
    cands = stage4_cands(c3, c3b, βC2)
    ds = combo_datasets(sys_obs(cands[1].sys))
    rows = fit_stage4(cands, ds)
    (id="s4", stage="4 Low resource + cisplatin coculture", title="Stage 4 · low resource + cisplatin coculture (1:1)",
     datasets=ds, u0=rows[1].u0, cands=rows, line="combination", kind=:s4)
end

# bootstrap with the rule replayed on every replicate (as in run_zoo_fits.jl)
function stage4_results(g; B)
    res = bootstrap(g; B=B)
    ks = [length(c.params) for c in g.cands]
    picks = size(res.bic, 1) > 0 ? [rule_pick(res.bic[b, :], ks) for b in 1:size(res.bic, 1)] : Int[]
    res = merge(res, (selected=1, pick=picks))
    (res=res, band=bands(g, res; tmax=XMAX, j=1))
end
