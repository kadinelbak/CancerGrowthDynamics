# Resistant cisplatin kill in coculture (PRD "Resistant-cell cisplatin kill in coculture", user 2026-10-23).
# Included by zoo_analyses.jl after the validation section (needs groups, gid, ck, save_ck, stale, g3, c3,
# βN2, βC2, evaluate, gfit, nominal_u0, nobs and the HTML helpers).
#
# A2780cis keeps its Stage 2b kill in coculture, β_C,co = λ_C·β_C. Three variants of the carried Stage 3a
# structure (same crowding and death terms, same residual weights) are compared:
#   zero       β_C,co = 0: resistant cells immune in coculture (the old scenario A, 2026-10-06 to 2026-10-22)
#   inherited  λ_C = 1: the Stage 2b kill unchanged
#   lambda     λ_C refitted
# 1. λ_C: estimate, profile-likelihood 95% interval (weighted SSE, χ²₁: n·log(SSE_w(λ)/SSE_w,min) ≤ 3.84),
#    likelihood-ratio test and ΔBIC against λ_C = 1 (the test run_zoo_fits.jl applies when it carries λ_C).
# 2. Held-out prediction with existing data only:
#    H1  leave one Ce1 seeding ratio out: refit on Ce0 + the other Ce1 ratios, predict it
#    H2  the low-resource + 1 µM coculture (Stage 4 data), predicted with no refit from each variant plus the
#        Stage 3b low-resource mechanisms (the "everything inherited" Stage 4 structure), nominal seeding
#    H3  the intermittent monoculture schedules: tests the inherited Stage 2b kill; the same for every variant
#    H4  crowding fitted to Ce0 only, Ce1 predicted (λ_C cannot be fitted without drug, so λ = 1 there)
#    Scores: weighted RMSE √(SSE_w/n) and RMSE in cells; "within 1 SE of the best" uses the SD of the paired
#    trajectory bootstrap of the weighted-RMSE difference.

include(joinpath(@__DIR__, "zoo_stage4.jl"))       # stage4_cands, combo_datasets (Stage 4 structure for H2)

const BC_GRID = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.0, 8.0]
const BC_VARIANTS = ("zero", "inherited", "lambda")
const BC_LABEL = Dict("zero" => "β_C,co = 0 (resistant cells immune in coculture; old scenario A)",
                      "inherited" => "λ_C = 1 (Stage 2b resistant kill inherited)", "lambda" => "λ_C refitted")

# the λ_C = 1 and λ_C-free candidates of the carried structure (stages.toml pairs every λ_C candidate with one)
function bc_pair()
    base = sort([s for s in c3.slots if s != :λC])
    find(withλ) = findfirst(r -> sort([s for s in r.slots if s != :λC]) == base && ((:λC in r.slots) == withλ), g3.cands)
    iB, iL = find(false), find(true)
    (iB === nothing || iL === nothing) && error("Stage 3a: the carried structure $(base) needs a λ_C = 1 and a λ_C-free candidate in stages.toml")
    g3.cands[iB], g3.cands[iL]
end
# a candidate with the resistant multiplier fixed at v (rebuilt from the λ_C = 1 candidate)
function with_lambda(cB, v)
    LN, LC = cB.sys.lines
    sys = build_system([LN, merge(LC, (s=merge(LC.s, (λ=Slot(0, v),)),))], cB.sys.envs, cB.sys.alpha)
    merge(cB, (sys=sys, model=make_system_model(sys), p0s=vcat([cB.params], cB.p0s), name="$(cB.name)_lambda$(v)"))
end
λidx(c) = findfirst(==(:λC), c.slots)
λval(c) = (i = λidx(c); i === nothing ? (hasproperty(c.sys.lines[2].s, :λ) ? c.sys.lines[2].s.λ.v : 1.0) : c.params[i])

function bc_profile(cB, cL)
    n = nobs(g3.datasets)
    fits = Vector{Any}(undef, length(BC_GRID))
    Threads.@threads :dynamic for j in eachindex(BC_GRID)
        fits[j] = try fit_cand(with_lambda(cB, BC_GRID[j]), g3.datasets) catch err; @warn "λ profile point failed" BC_GRID[j] err; nothing end
    end
    pts = [(λ=BC_GRID[j], sse=f.sse) for (j, f) in enumerate(fits) if f !== nothing && isfinite(f.sse)]
    push!(pts, (λ=λval(cL), sse=cL.sse)); sort!(pts; by=q -> q.λ)
    smin = minimum(q.sse for q in pts)
    dev = [(λ=q.λ, dev=n * log(q.sse / smin)) for q in pts]
    thr = 3.841458820694124
    λh = dev[argmin([q.dev for q in dev])].λ
    function cross(side)
        seq = side < 0 ? reverse([q for q in dev if q.λ <= λh]) : [q for q in dev if q.λ >= λh]
        for k in 2:length(seq)
            a, b = seq[k-1], seq[k]
            b.dev >= thr && return a.λ + (thr - a.dev) * (b.λ - a.λ) / (b.dev - a.dev)
        end
        side < 0 ? -Inf : Inf                                 # not reached inside the grid (or at the bound 0)
    end
    lo = cross(-1); hi = cross(1)
    stat = n * log(cB.sse / cL.sse)
    (points=dev, est=λval(cL), lo=lo == -Inf ? 0.0 : lo, lo_open=lo == -Inf, hi=hi, hi_open=hi == Inf, threshold=thr,
     lrt=stat, p=chisq1_sf(stat), dbic=cB.bic - cL.bic, n=n, fits=fits)
end
# P(χ²₁ > x) (as in run_zoo_fits.jl)
function chisq1_sf(x)
    x <= 0 && return 1.0
    z = sqrt(x / 2); t = 1 / (1 + 0.3275911z)
    t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * exp(-z^2)
end

# per-trajectory weighted SSE and point counts of a prediction
traj_wsse(ev, ds) = [sum(((p - y) / w)^2 for (p, y, w) in zip(P, d.y, dweights(d))) for (P, d) in zip(ev.predictions, ds)]
# score of every variant on one held-out set, with the paired trajectory bootstrap SD of the wRMSE difference to the best
function bc_score(evs::Dict, ds; B=2000, seed=20261023)
    nn = [length(d.x) for d in ds]; m = length(ds)
    S = Dict(v => traj_wsse(ev, ds) for (v, ev) in evs)
    w = Dict(v => sqrt(sum(S[v]) / sum(nn)) for v in keys(evs))
    best = argmin(w)
    rng = Xoshiro(seed); diffs = Dict(v => zeros(B) for v in keys(evs))
    for b in 1:B
        ix = rand(rng, 1:m, m); n = sum(nn[ix])
        wb = sqrt(sum(S[best][ix]) / n)
        for v in keys(evs)
            diffs[v][b] = sqrt(sum(S[v][ix]) / n) - wb
        end
    end
    Dict(v => (wrmse=w[v], rmse=evs[v].rmse, n=sum(nn), best=(v == best), diff=w[v] - w[best], se=std(diffs[v]),
               within=(v == best) || (w[v] - w[best] <= std(diffs[v]))) for v in keys(evs))
end

const BC_TAGS = dataset_tags("s3")
function bc_heldout(V)
    cB, cL = bc_pair()
    cZ = fit_cand(with_lambda(cB, 0.0), g3.datasets)
    C = Dict("zero" => cZ, "inherited" => cB, "lambda" => cL)
    D = g3.datasets
    # H1: leave one Ce1 ratio out
    ce1 = [k for (k, s) in enumerate(BC_TAGS) if string(get(s, "condition", "")) == "Ce1"]
    ratios = unique(string(BC_TAGS[k]["ratio"]) for k in ce1)
    h1 = Dict{String,Any}()
    for r in ratios
        held = [k for k in ce1 if string(BC_TAGS[k]["ratio"]) == r]; train = setdiff(eachindex(D), held)
        evs = Dict{String,Any}()
        jobs = collect(BC_VARIANTS); out = Vector{Any}(undef, length(jobs))
        Threads.@threads :dynamic for j in eachindex(jobs)
            c = C[jobs[j]]
            f = gfit(c.model, D[train], c.u0, c.params, c.bounds; starts=3, seed=7)
            out[j] = evaluate(c.model, D[held], nominal_u0(c, [D[k].state_index for k in held]), f.params)
        end
        for (j, v) in enumerate(jobs); evs[v] = out[j]; end
        h1[r] = (held=held, evs=evs)
    end
    held1 = reduce(vcat, [h1[r].held for r in ratios])
    ev1 = Dict(v => (predictions=reduce(vcat, [h1[r].evs[v].predictions for r in ratios]),
                     rmse=sqrt(sum(h1[r].evs[v].raw_sse for r in ratios) / sum(h1[r].evs[v].n_points for r in ratios))) for v in BC_VARIANTS)
    s1 = bc_score(ev1, D[held1])
    # H2: the combination culture, everything inherited (no refit), nominal seeding
    c3b = gid("s3b").cands[1]; rk = gid("s2-a2780cis").cands[1]; βC2v = rk.params[findfirst(==(:β), rk.slots)]
    ev2 = Dict{String,Any}(); ds2 = nothing
    for v in BC_VARIANTS
        c4 = only(c for c in stage4_cands(C[v], c3b, βC2v) if c.name == "s4_inherited")
        ds2 = combo_datasets(sys_obs(c4.sys))
        ev2[v] = evaluate(c4.model, ds2, c4.u0, Float64[])
    end
    s2 = bc_score(ev2, ds2)
    # H4: crowding fitted to Ce0 (the Validation page's held-out Ce1), resistant kill 0 or Stage 2b
    p4 = V.ce1[1].params; D1 = D[CE1]
    ev4 = Dict("zero" => evaluate(heldout_ce1_model(βN2, 0.0), D1, U0_CE1, p4), "inherited" => V.ce1[1].pred)
    s4 = bc_score(ev4, D1)
    (cands=C, h1=h1, ratios=ratios, s1=s1, ev2=ev2, s2=s2, ev4=ev4, s4=s4, βC2=βC2v)
end

if !haskey(ck, :betaC) || "s3" in stale || "s3b" in stale
    println("\n--- resistant kill in coculture: λ_C profile and held-out comparison ---"); flush(stdout)
    t0 = time()
    cB, cL = bc_pair()
    prof = bc_profile(cB, cL)
    ho = bc_heldout(ck[:validation])
    ck[:betaC] = (profile=prof, heldout=ho, carried=c3.name, pair=(cB.name, cL.name)); save_ck()
    @printf("  λ_C = %.3g (95%% profile %.3g–%s), LRT p = %.3g, ΔBIC(λ=1 − λ free) = %+.2f [%.0fs]\n", prof.est, prof.lo,
            prof.hi_open ? "∞" : @sprintf("%.3g", prof.hi), prof.p, prof.dbic, time() - t0); flush(stdout)
end

# ── page ──────────────────────────────────────────────────────────────────────
function profile_lambda_svg(P; W=560, H=260)
    pts = [q for q in P.points if isfinite(q.dev)]
    xmax = maximum(q.λ for q in pts); ymax = max(maximum(min(q.dev, 40.0) for q in pts), 6.0)
    X(x) = 50 + (W - 70) * x / xmax; Y(y) = H - 35 - (H - 60) * min(y, 40.0) / ymax
    path = join(["$(k == 1 ? "M" : "L")$(round(X(q.λ); digits=1)),$(round(Y(q.dev); digits=1))" for (k, q) in enumerate(pts)], " ")
    ticks = join(["<text x='$(X(x))' y='$(H - 15)' font-size='11' text-anchor='middle'>$(x)</text>" for x in (0, 1, 2, 3, 5, 8) if x <= xmax])
    """<svg viewBox='0 0 $W $H'><rect x='50' y='25' width='$(W - 70)' height='$(H - 60)' fill='#fff' stroke='#d3dde2'/>
<line x1='50' x2='$(W - 20)' y1='$(Y(P.threshold))' y2='$(Y(P.threshold))' stroke='#b45309' stroke-dasharray='4 3'/>
<line x1='$(X(1.0))' x2='$(X(1.0))' y1='25' y2='$(H - 35)' stroke='#5b6b73' stroke-dasharray='2 3'/>
<path d='$path' fill='none' stroke='$BLUE' stroke-width='2'/>$(join(["<circle cx='$(X(q.λ))' cy='$(Y(q.dev))' r='3' fill='$BLUE'/>" for q in pts]))
$ticks<text x='$(W / 2)' y='$(H - 1)' font-size='12' text-anchor='middle'>λ_C (coculture multiplier of the Stage 2b resistant kill)</text>
<text x='14' y='$(H / 2)' font-size='12' text-anchor='middle' transform='rotate(-90 14 $(H / 2))'>n·log(SSE_w/SSE_w,min)</text>
<text x='$(W - 24)' y='$(Y(P.threshold) - 4)' font-size='11' text-anchor='end' fill='#b45309'>95% threshold 3.84</text>
<text x='$(X(1.0) + 4)' y='38' font-size='11' fill='#5b6b73'>λ = 1</text></svg>"""
end

function write_resistant_kill()
    haskey(ck, :betaC) || return
    R = ck[:betaC]; P = R.profile; H = R.heldout
    f3(x) = isfinite(x) ? @sprintf("%.3g", x) : "—"
    hi = P.hi_open ? "> $(BC_GRID[end]) (not reached)" : f3(P.hi)
    lo = P.lo_open ? "0 (bound)" : f3(P.lo)
    βC = H.βC2
    carried = gid("s3").cands[1]
    carried_λ = λval(carried)
    carried_var = :λC in carried.slots ? "lambda" : "inherited"
    lrt = get(gid("s3"), :lrt, nothing)
    sets = [("H1 · leave one Ce1 ratio out", H.s1, "Ce0 + the other two Ce1 ratios"), ("H2 · low resource + 1 µM coculture (no refit)", H.s2, "nothing (Stage 3a + Stage 3b inherited)"),
            ("H4 · Ce1 from crowding fitted to Ce0", H.s4, "Ce0 only")]
    cell(S, v) = haskey(S, v) ? (S[v].best ? "<b>$(@sprintf("%.3f", S[v].wrmse))</b>" : @sprintf("%.3f", S[v].wrmse)) *
                                "<br><span class='muted'>$(fmt(S[v].rmse)) cells$(S[v].best ? "" : "; +" * @sprintf("%.3f", S[v].diff) * " ± " * @sprintf("%.3f", S[v].se))</span>" : "= λ_C = 1 here"
    rows = [trow([td(hesc(n)), td(hesc(tr)), tn(string(first(values(S)).n)), [td(cell(S, v)) for v in BC_VARIANTS]...]) for (n, S, tr) in sets]
    I = get(ck, :intermittent, nothing)
    if I !== nothing
        r = I["A2780cis"]
        push!(rows, trow([td("H3 · intermittent schedules, A2780cis monoculture (no refit)"), td("nothing (Stage 2b)"), tn(string(r.n)),
                          td("not applicable (monoculture)"), td("$(fmt(r.rmse)) cells (nRMSE $(@sprintf("%.2f", r.nrmse)))"), td("the same (λ_C acts in coculture only)")]))
    end
    # pass count: best or within 1 SE on each discriminating set; H3 counts for both variants that use the Stage 2b kill
    pass(v) = count(S -> haskey(S, v) ? S[v].within : S["inherited"].within, (H.s1, H.s2, H.s4)) + (I !== nothing && v != "zero" ? 1 : 0)
    passrows = [trow([td(hesc(BC_LABEL[v])), tn("$(pass(v)) of $(I === nothing ? 3 : 4)")]; cls=v == carried_var ? "win" : "") for v in BC_VARIANTS]
    pr = [trow([tn(f3(q.λ)), tn(@sprintf("%.2f", q.dev))]) for q in P.points]
    ntot = I === nothing ? 3 : 4
    good = [v for v in BC_VARIANTS if pass(v) >= 3]
    verdict = isempty(good) ?
        "<p class='note'><b>No variant passes $(min(3, ntot)) of $ntot.</b> The held-out sets disagree: $(join(["$(n): best $(hesc(BC_LABEL[only(v for v in keys(S) if S[v].best)]))" for (n, S, _) in sets], "; ")). The resistant kill in coculture is therefore not settled by the existing data. Read control results only where a schedule is robust across the bootstrap spread of λ_C·β_C, and treat both λ_C ≈ 0 and λ_C = 1 as plausible.</p>" :
        "<p class='note'>Passing at least 3 of $ntot: <b>$(join([hesc(BC_LABEL[v]) for v in good], "; "))</b>.$(carried_var in good ? " The carried model is among them." : " The carried model is not.")</p>"
    page("resistant_kill", "Resistant cisplatin kill in coculture", "resistant_kill", """
<div class='intro'><p>Does cisplatin still kill A2780cis cells when they grow with A2780? The coculture model now keeps the resistant line's Stage 2b kill and scales it by one coculture multiplier: \\(\\beta_{C,co} = \\lambda_C\\,\\beta_C\\), with \\(\\beta_C = $(fmt(βC))\\) from the A2780cis monoculture fit (Stage 2b). \\(\\lambda_C = 1\\) means the monoculture kill carries over unchanged; \\(\\lambda_C = 0\\) means resistant cells are immune in coculture (the old scenario A). The old fixed-zero candidates no longer compete in Stage 3a; \\(\\beta_{C,co} = 0\\) is kept here only as a comparison (decision: PRD of 2026-10-23, which reverses the 2026-10-06 choice).</p>
<p>Carried Stage 3a model: <b>$(hesc(carried.label))</b>, so the resistant kill in coculture is \\(\\beta_{C,co} = $(f3(carried_λ)) \\times $(fmt(βC)) = $(fmt(carried_λ * βC))\\)$(lrt === nothing ? "" : lrt.keep ? " (λ_C free: picked by the selection rule and the likelihood-ratio test, p = $(f3(lrt.p)))" : " (the rule picked λ_C free, but the likelihood-ratio test did not support it, p = $(f3(lrt.p)), so λ_C = 1 is carried)").</p></div>
<h2>How strongly the coculture data pin λ_C down $(tip("Profile likelihood: λ_C is fixed at each value and every other Stage 3a parameter is refitted (weighted, as the carried fits). The curve is n·log(SSE_w(λ)/SSE_w,min); values below the dashed line (3.84, χ² with 1 degree of freedom) form the 95% interval.", "Narrower is better: the coculture data then determine λ_C."))</h2>
<div class='grid'><div class='panel'>$(profile_lambda_svg(P))</div><div class='panel'>
<table><tr><th>λ_C estimate</th><td class='num'>$(f3(P.est))</td></tr><tr><th>95% profile interval</th><td class='num'>$lo to $hi</td></tr>
<tr><th>β_C,co = λ_C·β_C</th><td class='num'>$(f3(P.est * βC)) (interval $(P.lo_open ? "0" : f3(P.lo * βC)) to $(P.hi_open ? "—" : f3(P.hi * βC)))</td></tr>
<tr><th>Likelihood-ratio test vs λ_C = 1</th><td class='num'>LR = $(@sprintf("%.2f", P.lrt)), p = $(f3(P.p))</td></tr>
<tr><th>ΔBIC (λ_C = 1 minus λ_C free)</th><td class='num'>$(@sprintf("%+.2f", P.dbic))</td></tr>
<tr><th>Structure compared</th><td>$(hesc(R.pair[1])) vs $(hesc(R.pair[2]))</td></tr></table>
<p class='muted'>Keep λ_C free only when ΔBIC &gt; 3 and p &lt; 0.05 (stages.toml <code>lrt_p</code>). The tested doses (0.62–1.47 µM) sit well below the A2780cis IC50 of 8.34 µM, so a wide interval is expected; when it is wide, λ_C = 1 is the default.</p></div></div>
<details><summary>Profile points</summary>$(table([("λ_C", "The fixed multiplier."), ("n·log(SSE ratio)", "Weighted fit loss relative to the best value; below 3.84 is inside the 95% interval.")], pr))</details>
<h2>Which resistant kill predicts unseen data best $(tip("Each variant has the same crowding and death terms and residual weights as the carried model and differs only in the resistant kill. The held-out data were not used to fit the prediction being scored.", :none))</h2>
$(table([("Held-out set", "What is predicted."), ("Fitted to", "The data the prediction was fitted to."), ("Points", "Measurements scored."),
         ("β_C,co = 0", "Weighted RMSE √(SSE_w/n) (bold = best), the RMSE in cells, and the difference to the best ± its paired-bootstrap SD over held-out trajectories."),
         ("λ_C = 1", "As the first score column, for the inherited Stage 2b resistant kill."), ("λ_C refitted", "As the first score column, with λ_C refitted on the training data.")], rows))
<p>Best or within one standard error of the best:</p>
$verdict
$(table([("Variant", "The resistant kill in coculture."), ("Passes", "Held-out sets on which the variant is best or within one SE of the best (H3 counts for the variants that use the Stage 2b kill). The PRD asks the carried model to pass at least 3 of 4.")], passrows))
<p class='muted'>H1 refits on 2 of the 3 Ce1 ratios; H2 uses data that Stage 4 is fitted to (here only Stages 3a and 3b are inherited, with no refit, so it is held out for this question); H3 does not separate the variants because λ_C acts only in coculture; H4 cannot fit λ_C (no drug in Ce0), so λ_C refitted equals λ_C = 1 there. Ce1 informed earlier modelling choices, so H4 is not fully independent. Interval coverage is not scored: the bootstrap bands of the three variants are not refitted here. Optimal control with this model, including how robust each schedule is across the β_C spread: <a href='../../control/report.html'>Optimal control</a>.</p>""")
end
