# Held-out treated-coculture test (validation.toml [heldout]) and the fitting / scoring helpers it shares
# with zoo_analyses.jl. Needs `groups`, `gid` and `GPE` (= GrowthParameterEstimation) defined first.
# Used by zoo_analyses.jl (Validation page) and run_heldout.jl (any zoo run, e.g. a sensitivity variant);
# both write the same outputs/…/heldout_ce1.csv through write_heldout_csv.
#
# Every refit here uses the run's residual weights, like the carried fits: datasets that carry `w`
# (stages.toml residual_weighting) go to the package through gpe_datasets (one dataset per point with its
# residual_scale), so the optimiser minimises SSE_w = Σ((ŷ − y)/w)², the criterion of the carried fits.
# Scores are reported as RMSE in cells (raw, for interpretation) and as weighted RMSE √(SSE_w/n).
# Starting values (u0) are fixed nuisance values; they are never refitted here. No wall-clock limits:
# every optimiser call stops on its iteration cap, so the results do not depend on the machine.

using Random: Xoshiro
using Statistics

const FIT_KW = (optimizer=:nelder_mead, reltol=1e-8, abstol=1e-8, initial_time=0.0)
const FIT_MAXITERS = 2000
unweighted(ds) = [(x=d.x, y=d.y, state_index=d.state_index) for d in ds]
nobs(ds) = sum(length(d.x) for d in ds)
# predictions per dataset, the raw SSE / RMSE in cells and the weighted SSE / RMSE (the fits' criterion)
function evaluate(model, ds, u0, p)
    ev = GPE.evaluate_joint_fit(model, unweighted(ds), u0, p; initial_time=0.0, reltol=1e-8, abstol=1e-8)
    wsse = sum(sum(((pr - y) / w)^2 for (pr, y, w) in zip(P, d.y, dweights(d))) for (d, P) in zip(ds, ev.predictions))
    (predictions=ev.predictions, raw_sse=ev.raw_sse, n_points=ev.n_points, rmse=ev.rmse, wsse=wsse, wrmse=sqrt(wsse / ev.n_points))
end
# weighted multistart refit (the supplied p0 is the first start); .sse of the result is SSE_w
function gfit(model, ds, u0, p0, bounds; starts=3, seed=1, maxiters=FIT_MAXITERS)
    st = GPE.generate_multistarts(p0, bounds; n_starts=starts, seed=seed)
    GPE.run_joint_multistart(model, gpe_datasets(collect(ds)), u0, st; bounds=bounds, maxiters=maxiters, FIT_KW...).fit
end

# Starting state with the nominal seeding for the given observed states. A candidate whose fit estimated
# some starting values carries x0est (indices into u0) and x0nom (nominal values: one per observed state
# in u0 order, or one per x0est entry); without x0nom its u0 is returned unchanged.
function nominal_u0(c, states)
    u = copy(c.u0)
    hasproperty(c, :x0nom) || return u
    nom = c.x0nom; est = hasproperty(c, :x0est) ? c.x0est : Int[]
    nob = hasproperty(c, :sys) ? length(sys_obs(c.sys)) : -1
    for s in states
        if length(nom) == nob
            u[s] = nom[s]
        elseif length(nom) == length(est) && (j = findfirst(==(s), est)) !== nothing
            u[s] = nom[j]
        end
    end
    u
end
estimated_starts(c) = hasproperty(c, :x0est) ? Set{Int}(c.x0est) : Set{Int}()

# validation.toml: which data are held out and how they are predicted
const VAL = TOML.parsefile(joinpath(@__DIR__, "validation.toml"))

# A held-out condition (validation.toml [heldout]): the group's crowding is refitted (weighted) on the
# training data only (the untreated Ce0 cocultures), the kill strengths come from the monoculture drug
# fits (Stage 3a's own were fitted with Ce1, so they are not used), and the held-out data are predicted
# from the nominal seeding (a starting value the full fit estimated from Ce1 is not used).
const HO = VAL["heldout"]
g3 = gid(HO["group"]); c3 = g3.cands[1]
βof(c) = c.params[findfirst(==(:β), c.slots)]
βN2 = βof(gid("$(HO["drug_from"])-a2780").cands[1]); βC2 = βof(gid("$(HO["drug_from"])-a2780cis").cands[1])
function heldout_ce1_sys(βN, βC)
    LN, LC = c3.sys.lines
    # λ_C = 1: the resistant kill is exactly the given βC (Stage 2b), whatever multiplier Stage 3a refitted
    LN2 = merge(LN, (s=merge(LN.s, (β=Slot(0, βN),)),)); LC2 = merge(LC, (s=merge(LC.s, (β=Slot(0, βC), λ=Slot(0, 1.0))),))
    build_system([LN2, LC2], c3.sys.envs, c3.sys.alpha)
end
heldout_ce1_model(βN, βC) = make_system_model(heldout_ce1_sys(βN, βC))
# largest parameter index any slot of a system reads (lines, their extra modules, crowding)
max_free_index(sys) = maximum(vcat([q.i for L in sys.lines for M in modules(L) for q in values(M.s) if q isa Slot],
                                   [q.i for row in sys.alpha for q in row if q isa Slot]); init=0)
tag_match(s, sel) = all(string(get(s, k, "")) == string(v) for (k, v) in sel)
const HO_TAGS = dataset_tags(HO["group"])
const CE0 = [k for (k, s) in enumerate(HO_TAGS) if tag_match(s, HO["train"])]
const CE1 = [k for (k, s) in enumerate(HO_TAGS) if tag_match(s, HO["predict"])]
ds_label(s) = join(filter(!isempty, [string(get(s, "condition", "")), string(get(s, "ratio", "")), string(get(s, "line", ""))]), " ")
# the listed crowding parameters the carried model has (e.g. a model without sensitive death has no d_N);
# they must lead its parameter vector, so the kill slots after them can be fixed to the monoculture values
const HO_NAMES = [n for n in HO["refit"] if n in c3.pnames]
const HO_RI = [findfirst(==(n), c3.pnames) for n in HO_NAMES]
isempty(HO_RI) && error("$(c3.name) has none of the validation.toml heldout.refit parameters $(HO["refit"])")
HO_RI == collect(1:length(HO_RI)) || error("validation.toml heldout.refit must be the first parameters of $(c3.name): $(c3.pnames)")
max_free_index(heldout_ce1_sys(βN2, βC2)) <= length(HO_RI) ||
    error("$(c3.name) has free parameters other than heldout.refit $(HO_NAMES) and the kills: $(c3.pnames); add them to validation.toml heldout.refit")
# Ce1 starting values: nominal seeding; the Ce1 series whose start the full fit had estimated are reported
const CE1_STATES = [g3.datasets[k].state_index for k in CE1]
const U0_CE1 = nominal_u0(c3, CE1_STATES)
const CE1_EST = [ds_label(HO_TAGS[k]) for k in CE1 if g3.datasets[k].state_index in estimated_starts(c3)]
const CE1_U0_NOTE = hasproperty(c3, :x0nom) ?
    "nominal seeding for every Ce1 series" * (isempty(CE1_EST) ? " (the full fit estimated none of them)" : " (the full fit had estimated the start of: $(join(CE1_EST, ", ")))") :
    "the carried model's starting values (it estimated none)"

# Paired bootstrap of RMSE(prediction) − RMSE(reference) over the held-out trajectories: trajectories are
# resampled with replacement (predictions fixed, no refit) and both RMSEs are pooled over the same draw.
const HO_BOOT = 10_000
function paired_rmse_boot(pred, ref, ds; B=HO_BOOT, seed=20261006)
    sp = [sum(abs2, P .- d.y) for (P, d) in zip(pred, ds)]; sr = [sum(abs2, P .- d.y) for (P, d) in zip(ref, ds)]
    nn = [length(d.y) for d in ds]; m = length(ds); rng = Xoshiro(seed)
    diffs = Vector{Float64}(undef, B)
    for b in 1:B
        ix = rand(rng, 1:m, m); n = sum(nn[ix])
        diffs[b] = sqrt(sum(sp[ix]) / n) - sqrt(sum(sr[ix]) / n)
    end
    (diff=sqrt(sum(sp) / sum(nn)) - sqrt(sum(sr) / sum(nn)), lo=quantile(diffs, 0.025), hi=quantile(diffs, 0.975),
     share_below0=mean(diffs .< 0), B=B, n_traj=m)
end

function heldout_ce1()
    m = heldout_ce1_model(βN2, βC2)
    D = g3.datasets                       # carry w when the run is weighted
    f = gfit(m, D[CE0], U0_CE1, c3.params[HO_RI], c3.bounds[HO_RI]; starts=4)
    pred = evaluate(m, D[CE1], U0_CE1, f.params)
    ce1 = [(variant="cisplatin from Stage 2 for both lines", params=f.params, fit_ce0=evaluate(m, D[CE0], U0_CE1, f.params), pred=pred)]
    full = evaluate(c3.model, D[CE1], c3.u0, c3.params)          # in-sample: the full fit with its own starts
    # reference with no drug model: the same Ce0-fitted growth and crowding, no drug term
    HO["baseline"] == "logistic" || error("validation.toml heldout.baseline: only \"logistic\" is implemented")
    base = evaluate(heldout_ce1_model(0.0, 0.0), D[CE1], U0_CE1, f.params)
    boot = paired_rmse_boot(pred.predictions, base.predictions, D[CE1])
    (ce1=ce1, ce1_full=full, ce1_baseline_rmse=base.rmse, ce1_baseline="logistic", ce1_baseline_ev=base, ce1_boot=boot,
     ce1_est_starts=CE1_EST, ce1_u0_note=CE1_U0_NOTE)
end

# outputs/…/heldout_ce1.csv, written by zoo_analyses.jl (from its checkpoint) and run_heldout.jl alike
function write_heldout_csv(r, dir)
    D1 = g3.datasets[CE1]; m1 = mean(reduce(vcat, [d.y for d in D1])); n1 = nobs(D1)
    p = r.ce1[1]; b = get(r, :ce1_boot, nothing); bev = get(r, :ce1_baseline_ev, nothing)
    note = get(r, :ce1_u0_note, "")
    row(name, rmse, wrmse, params, lo, hi, nt) = (row=name, rmse=rmse, nrmse=rmse / m1, wrmse=wrmse, n=n1, params=params, diff_lo=lo, diff_hi=hi, note=nt)
    rows = [row("prediction (crowding from Ce0, kills from Stage 2b)", p.pred.rmse, get(p.pred, :wrmse, NaN),
                join(["$(n)=$(v)" for (n, v) in zip(c3.pnames[HO_RI], p.params)], "; "), NaN, NaN, "Ce1 starts: " * note),
            row("in-sample (Stage 3a fitted with Ce1)", r.ce1_full.rmse, get(r.ce1_full, :wrmse, NaN), "", NaN, NaN, ""),
            row("untreated reference (same crowding, no drug)", r.ce1_baseline_rmse, bev === nothing ? NaN : bev.wrmse, "", NaN, NaN, "")]
    b === nothing || push!(rows, row("prediction minus untreated reference (RMSE difference)", b.diff, NaN, "", b.lo, b.hi,
                                     "95% paired bootstrap over the $(b.n_traj) Ce1 trajectories, B = $(b.B); share of draws < 0: $(round(b.share_below0; digits=3))"))
    CSV.write(joinpath(dir, "heldout_ce1.csv"), DataFrame(rows))
    println("Wrote $(joinpath(dir, "heldout_ce1.csv"))")
end
