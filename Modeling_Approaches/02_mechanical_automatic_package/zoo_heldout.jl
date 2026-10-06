# Held-out treated-coculture test (validation.toml [heldout]) and the fitting helpers it shares with
# zoo_analyses.jl. Needs `groups`, `gid` and `GPE` (= GrowthParameterEstimation) defined first.
# Used by zoo_analyses.jl (Validation page) and run_heldout.jl (any zoo run, e.g. a sensitivity variant).
# The held-out fit is unweighted whatever the run's residual_weighting, so RMSEs compare across runs.

const FIT_KW = (optimizer=:nelder_mead, maxiters=2000, maxtime=180.0, reltol=1e-8, abstol=1e-8, initial_time=0.0)
evaluate(model, ds, u0, p) = GPE.evaluate_joint_fit(model, ds, u0, p; initial_time=0.0, reltol=1e-8, abstol=1e-8)
function gfit(model, ds, u0, p0, bounds; starts=3, seed=1)
    st = GPE.generate_multistarts(p0, bounds; n_starts=starts, seed=seed)
    GPE.run_joint_multistart(model, ds, u0, st; bounds=bounds, FIT_KW...).fit
end
nobs(ds) = sum(length(d.x) for d in ds)
unweighted(ds) = [(x=d.x, y=d.y, state_index=d.state_index) for d in ds]

# validation.toml: which data are held out and how they are predicted
const VAL = TOML.parsefile(joinpath(@__DIR__, "validation.toml"))

# A held-out condition (validation.toml [heldout]): the group's crowding is refitted on the training
# data only (the untreated Ce0 cocultures), the kill strengths come from the monoculture drug fits
# (Stage 3a's own were fitted with Ce1, so they are not used), and the held-out data are predicted.
const HO = VAL["heldout"]
g3 = gid(HO["group"]); c3 = g3.cands[1]
βof(c) = c.params[findfirst(==(:β), c.slots)]
βN2 = βof(gid("$(HO["drug_from"])-a2780").cands[1]); βC2 = βof(gid("$(HO["drug_from"])-a2780cis").cands[1])
function heldout_ce1_model(βN, βC)
    LN, LC = c3.sys.lines
    LN2 = merge(LN, (s=merge(LN.s, (β=Slot(0, βN),)),)); LC2 = merge(LC, (s=merge(LC.s, (β=Slot(0, βC),)),))
    sys = build_system([LN2, LC2], c3.sys.envs, c3.sys.alpha)
    make_system_model(sys)
end
tag_match(s, sel) = all(string(get(s, k, "")) == string(v) for (k, v) in sel)
const HO_TAGS = dataset_tags(HO["group"])
const CE0 = [k for (k, s) in enumerate(HO_TAGS) if tag_match(s, HO["train"])]
const CE1 = [k for (k, s) in enumerate(HO_TAGS) if tag_match(s, HO["predict"])]
# the refitted parameters must lead the carried model's parameter vector (the kill slots are then fixed)
# the listed crowding parameters the carried model has (e.g. a model without sensitive death has no d_N);
# they must lead its parameter vector, so the kill slots after them can be fixed to the monoculture values
const HO_NAMES = [n for n in HO["refit"] if n in c3.pnames]
const HO_RI = [findfirst(==(n), c3.pnames) for n in HO_NAMES]
isempty(HO_RI) && error("$(c3.name) has none of the validation.toml heldout.refit parameters $(HO["refit"])")
HO_RI == collect(1:length(HO_RI)) || error("validation.toml heldout.refit must be the first parameters of $(c3.name): $(c3.pnames)")
function heldout_ce1()
    m = heldout_ce1_model(βN2, βC2)
    D = unweighted(g3.datasets)
    f = gfit(m, D[CE0], c3.u0, c3.params[HO_RI], c3.bounds[HO_RI]; starts=4)
    ce1 = [(variant="cisplatin from Stage 2 for both lines", params=f.params, fit_ce0=evaluate(m, D[CE0], c3.u0, f.params),
            pred=evaluate(m, D[CE1], c3.u0, f.params))]
    full = evaluate(c3.model, D[CE1], c3.u0, c3.params)
    # reference with no drug model: the same Ce0-fitted growth and crowding, no drug term
    HO["baseline"] == "logistic" || error("validation.toml heldout.baseline: only \"logistic\" is implemented")
    base = evaluate(heldout_ce1_model(0.0, 0.0), D[CE1], c3.u0, f.params)
    (ce1=ce1, ce1_full=full, ce1_baseline_rmse=base.rmse, ce1_baseline="logistic")
end

