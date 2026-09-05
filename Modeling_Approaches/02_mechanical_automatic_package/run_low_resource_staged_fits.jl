#!/usr/bin/env julia
"""Fit the available low-resource trajectories with GrowthParameterEstimation.

Inputs are the reproducible workbook export made by `extract_low_resource_workbooks.py`.
The report intentionally distinguishes a shared Run-1/Run-2 *joint* shape fit
(normalised at each run's observed day zero) from a pooled mean trajectory.
"""

using CSV, DataFrames, Dates
using GrowthParameterEstimation
const GPE = GrowthParameterEstimation

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(@__DIR__, "outputs", "low_resource_staged_fits")
const INPUT = joinpath(@__DIR__, "outputs", "low_resource_workbook_trajectories.csv")
mkpath(OUT)

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
slug(s) = replace(lowercase(string(s)), r"[^a-z0-9]+" => "-")

function logistic!(du, u, p, t)
    r, K = p
    du[1] = r*u[1]*(1-u[1]/max(K, 1e-8))
end
function loss_logistic!(du, u, p, t)
    r, K, d = p
    du[1] = r*u[1]*(1-u[1]/max(K, 1e-8)) - d*u[1]
end
function delay_logistic!(du, u, p, t)
    r, K, lag = p
    du[1] = (t >= lag ? r : 0.0)*u[1]*(1-u[1]/max(K, 1e-8))
end

function candidate_specs(scale)
    cap = max(2.0, scale * 1.25)
    return Dict(
        "logistic" => (model=logistic!, p0=[0.35, cap], bounds=[(1e-6, 4.0), (1e-3, max(10.0, scale*20))]),
        "logistic_plus_loss" => (model=loss_logistic!, p0=[0.5, cap, 0.08], bounds=[(1e-6, 4.0), (1e-3, max(10.0, scale*20)), (0.0, 3.0)]),
        "logistic_with_lag" => (model=delay_logistic!, p0=[0.5, cap, 0.5], bounds=[(1e-6, 4.0), (1e-3, max(10.0, scale*20)), (0.0, 8.0)]),
    )
end

function fit_single(x, y)
    specs = candidate_specs(maximum(y))
    rows = NamedTuple[]; fits = Dict{String,Any}()
    for (name, spec) in specs
        # The package's joint API supports its bounded derivative-free optimiser.
        # This keeps the lag model in the model zoo even when its discontinuity
        # makes the default BFGS line search non-finite.
        try
            result = GPE.run_joint_fit(spec.model, [(x=x, y=y, state_index=1)], [y[1]], spec.p0;
                bounds=spec.bounds, maxiters=2_000, optimizer=:nelder_mead, reltol=1e-7, abstol=1e-7)
            fits[name] = result
            push!(rows, (model=name, bic=result.bic, sse=result.sse, params=join(round.(result.params; sigdigits=5), ";")))
        catch err
            # Preserve failed candidates in the BIC table rather than silently
            # removing a model from consideration.
            push!(rows, (model=name, bic=Inf, sse=Inf, params="fit failed: $(typeof(err))"))
        end
    end
    sort!(rows; by = row -> row.bic)
    return rows, fits
end

function fit_joint(runs)
    # Each run is scaled to its own day-zero value: shared r/K/loss describe shape,
    # while plotting restores original cell-count units for each run.
    specs = candidate_specs(maximum(maximum(run.y ./ run.y[1]) for run in runs))
    data = [(x=run.x, y=run.y ./ run.y[1], state_index=1) for run in runs]
    rows = NamedTuple[]; fits = Dict{String,Any}()
    for (name, spec) in specs
        result = GPE.run_joint_fit(spec.model, data, [1.0], spec.p0; bounds=spec.bounds, maxiters=2_000, optimizer=:nelder_mead, reltol=1e-7, abstol=1e-7)
        fits[name] = result
        push!(rows, (model=name, bic=result.bic, sse=result.sse, params=join(round.(result.params; sigdigits=5), ";")))
    end
    sort!(rows; by = row -> row.bic)
    return rows, fits
end

function svg_panel(title, series, predictions; width=720, height=300)
    xs = vcat([s.x for s in series]...); ys = vcat([s.y for s in series]...)
    ymax = max(maximum(ys), 1.0) * 1.10; xmin, xmax = minimum(xs), maximum(xs)
    sx(x) = 55 + (width-75)*(x-xmin)/max(xmax-xmin, 1e-8)
    sy(y) = height-35 - (height-60)*y/ymax
    colors = ["#0b6e99", "#d95f02", "#3b8b4f", "#9b59b6"]
    parts = ["<svg viewBox='0 0 $width $height' role='img'><rect width='100%' height='100%' fill='white'/><text x='55' y='18' font-size='14' font-weight='bold'>$(esc(title))</text><line x1='55' y1='$(height-35)' x2='$(width-20)' y2='$(height-35)' stroke='#555'/><line x1='55' y1='30' x2='55' y2='$(height-35)' stroke='#555'/><text x='8' y='36' font-size='10'>$(round(ymax; sigdigits=3))</text><text x='24' y='$(height-38)' font-size='10'>0</text>"]
    for (i, s) in enumerate(series)
        color = colors[mod1(i, length(colors))]
        pred = predictions[i]
        line = join(["$(round(sx(x); digits=1)),$(round(sy(y); digits=1))" for (x,y) in zip(s.x,pred)], " ")
        dots = join(["<circle cx='$(round(sx(x);digits=1))' cy='$(round(sy(y);digits=1))' r='3' fill='$color'/>" for (x,y) in zip(s.x,s.y)], "")
        push!(parts, "<polyline points='$line' fill='none' stroke='$color' stroke-width='2'/>$dots<text x='$(70+(i-1)*150)' y='$(height-8)' font-size='10' fill='$color'>● $(esc(s.name))</text>")
    end
    push!(parts, "<text x='$(width-55)' y='$(height-18)' font-size='10'>day</text></svg>")
    return join(parts)
end

function bic_svg(rows; width=560, height=170)
    base = minimum(r.bic for r in rows); span = max(maximum(r.bic for r in rows)-base, 1.0)
    pieces = ["<svg viewBox='0 0 $width $height'><text x='10' y='16' font-size='13' font-weight='bold'>BIC (lower is preferred)</text>"]
    for (i, row) in enumerate(rows)
        y = 28 + (i-1)*42; w = 60 + 430*(row.bic-base)/span
        push!(pieces, "<text x='10' y='$(y+14)' font-size='11'>$(esc(row.model))</text><rect x='170' y='$y' width='$w' height='23' fill='#3478a8'/><text x='$(175+w)' y='$(y+16)' font-size='11'>$(round(row.bic; digits=1))</text>")
    end
    push!(pieces, "</svg>"); return join(pieces)
end

raw = CSV.read(INPUT, DataFrame)
means = filter(:series => ==("Mean Cells"), raw)
mono = filter(row -> !startswith(lowercase(row.workbook), "ce"), means)
ce = filter(row -> startswith(lowercase(row.workbook), "ce"), means)

all_rows = NamedTuple[]; panels = String[]

# Stage LR-1: every available single-population workbook trajectory; run labels are
# derived from sheet names.  Shared run fits are produced only where both runs occur.
groups = groupby(mono, [:workbook, :sheet])
run_buckets = Dict{String,Vector{NamedTuple}}()
for g in groups
    x = Float64.(g.day); y = Float64.(g.value)
    order = sortperm(x); x, y = x[order], y[order]
    name = "$(g.workbook[1]) — $(g.sheet[1])"
    run = startswith(g.sheet[1], "1 -") ? "Run 1" : startswith(g.sheet[1], "2 -") ? "Run 2" : "single run"
    key = replace(name, r"^.*? — [12] - " => "")
    ranked, fits = fit_single(x, y); best = ranked[1]; pred = fits[best.model].predictions[1]
    append!(all_rows, [(stage="LR-1 monoculture/lineage", condition=name, fit_scope=run, model=r.model, bic=r.bic, sse=r.sse, params=r.params, inherited="none; direct low-resource fit") for r in ranked])
    push!(get!(run_buckets, key, NamedTuple[]), (name=run, x=x, y=y, best=best.model, fit=fits[best.model]))
    push!(panels, "<section><h3>$(esc(name)) — $(esc(run))</h3>$(svg_panel("Observed points and winning fit: $(best.model)", [(name=run,x=x,y=y)], [pred]))$(bic_svg(ranked))</section>")
end

for (key, runs) in run_buckets
    length(runs) == 2 || continue
    ranked, fits = fit_joint(runs); best = ranked[1]
    pred = [[fits[best.model].predictions[i][j] * runs[i].y[1] for j in eachindex(runs[i].x)] for i in eachindex(runs)]
    append!(all_rows, [(stage="LR-1 monoculture/lineage", condition=key, fit_scope="Run 1 + Run 2 (joint normalized shape)", model=r.model, bic=r.bic, sse=r.sse, params=r.params, inherited="shared r/K/loss across runs; run-specific observed N0") for r in ranked])
    push!(panels, "<section><h3>$(esc(key)) — joint Run 1 + Run 2</h3>$(svg_panel("Joint fitted shape restored to each run's observed N0", [(name=r.name,x=r.x,y=r.y) for r in runs], pred))$(bic_svg(ranked))</section>")
end

# Stage LR-2/LR-3: Ce0 provides untreated coculture baselines.  Ce1 is fitted with
# those baseline growth/competition parameters inherited; only treatment losses are
# fitted.  A single treatment level cannot identify Hill/IC50, so those are excluded.
function null2!(du,u,p,t)
    rS,KS,rR,KR = p; du[1]=rS*u[1]*(1-u[1]/max(KS,1e-8)); du[2]=rR*u[2]*(1-u[2]/max(KR,1e-8))
end
function comp2!(du,u,p,t)
    rS,KS,aSR,rR,KR,aRS=p; du[1]=rS*u[1]*(1-(u[1]+aSR*u[2])/max(KS,1e-8)); du[2]=rR*u[2]*(1-(u[2]+aRS*u[1])/max(KR,1e-8))
end
function ce_pair(g)
    ss = filter(row -> occursin("Sensitive", row.sheet), g); rr = filter(row -> occursin("Resistant", row.sheet), g)
    x=Float64.(ss.day); ix=sortperm(x); return x[ix], Float64.(ss.value)[ix], Float64.(rr.value)[ix]
end
ce0 = Dict{String,Any}()
for g in groupby(filter(row -> startswith(lowercase(row.workbook), "ce0"), ce), :workbook)
    x,S,R=ce_pair(g); scale=max(maximum(S),maximum(R)); data=[(x=x,y=S,state_index=1),(x=x,y=R,state_index=2)]
    specs=Dict("independent_logistic"=>(null2!,[0.4,scale,0.4,scale],[(1e-6,4.0),(1e-3,scale*20),(1e-6,4.0),(1e-3,scale*20)]), "competition_logistic"=>(comp2!,[0.4,scale,1.0,0.4,scale,1.0],[(1e-6,4.0),(1e-3,scale*20),(0.0,5.0),(1e-6,4.0),(1e-3,scale*20),(0.0,5.0)]))
    rows=NamedTuple[]; fits=Dict{String,Any}()
    for (model,(fun,p0,bounds)) in specs
        f=GPE.run_joint_fit(fun,data,[S[1],R[1]],p0;bounds=bounds,maxiters=3000,optimizer=:nelder_mead,reltol=1e-7,abstol=1e-7); fits[model]=f; push!(rows,(model=model,bic=f.bic,sse=f.sse,params=join(round.(f.params;sigdigits=5),";")))
    end
    sort!(rows; by = row -> row.bic); best=rows[1]; ce0[replace(g.workbook[1],"ce0_"=>"")]=(best=best, fit=fits[best.model], model=best.model)
    append!(all_rows,[(stage="LR-2 untreated coculture",condition=g.workbook[1],fit_scope="available ratio",model=r.model,bic=r.bic,sse=r.sse,params=r.params,inherited="lineage starting counts observed") for r in rows])
    pred=fits[best.model].predictions; push!(panels,"<section><h3>$(esc(g.workbook[1])) — untreated coculture</h3>$(svg_panel("Observed and winning inherited-baseline candidate",[(name="Sensitive",x=x,y=S),(name="Resistant",x=x,y=R)],pred))$(bic_svg(rows))</section>")
end

for g in groupby(filter(row -> startswith(lowercase(row.workbook), "ce1"), ce), :workbook)
    ratio=replace(g.workbook[1],"ce1_"=>""); haskey(ce0,ratio) || continue
    x,S,R=ce_pair(g); base=ce0[ratio].fit.params; data=[(x=x,y=S,state_index=1),(x=x,y=R,state_index=2)]
    inherited_fun = if ce0[ratio].model == "competition_logistic"
        (du,u,p,t)->begin rS,KS,aSR,rR,KR,aRS=base; dS,dR=p; du[1]=rS*u[1]*(1-(u[1]+aSR*u[2])/max(KS,1e-8))-dS*u[1]; du[2]=rR*u[2]*(1-(u[2]+aRS*u[1])/max(KR,1e-8))-dR*u[2] end
    else
        (du,u,p,t)->begin rS,KS,rR,KR=base; dS,dR=p; du[1]=rS*u[1]*(1-u[1]/max(KS,1e-8))-dS*u[1]; du[2]=rR*u[2]*(1-u[2]/max(KR,1e-8))-dR*u[2] end
    end
    f=GPE.run_joint_fit(inherited_fun,data,[S[1],R[1]],[0.1,0.1];bounds=[(0.0,3.0),(0.0,3.0)],maxiters=3000,optimizer=:nelder_mead,reltol=1e-7,abstol=1e-7)
    rows=[(model="inherited_$(ce0[ratio].model)_plus_lineage_treatment_loss",bic=f.bic,sse=f.sse,params=join(round.(f.params;sigdigits=5),";"))]
    append!(all_rows,[(stage="LR-3 treated coculture",condition=g.workbook[1],fit_scope="available ratio",model=r.model,bic=r.bic,sse=r.sse,params=r.params,inherited="$(ce0[ratio].model) parameters from ce0_$(ratio); fitted dS,dR only") for r in rows])
    push!(panels,"<section><h3>$(esc(g.workbook[1])) — treated coculture</h3>$(svg_panel("Observed and inherited untreated-coculture fit",[(name="Sensitive",x=x,y=S),(name="Resistant",x=x,y=R)],f.predictions))$(bic_svg(rows))</section>")
end

results=DataFrame(all_rows); CSV.write(joinpath(OUT,"bic_model_ranking.csv"),results)
summary_rows = join(["<tr><td>$(esc(r.stage))</td><td>$(esc(r.condition))</td><td>$(esc(r.fit_scope))</td><td>$(esc(r.model))</td><td>$(round(r.bic;digits=2))</td><td>$(esc(r.inherited))</td></tr>" for r in eachrow(results) if r.bic == minimum(results.bic[(results.stage .== r.stage) .& (results.condition .== r.condition) .& (results.fit_scope .== r.fit_scope)])], "\n")
html = """<!doctype html><html><head><meta charset='utf-8'><title>Low-resource staged parameter estimation</title><style>body{font:14px system-ui;margin:32px;color:#17212b;max-width:1250px}h1{color:#124b6e}section{border-top:1px solid #ccd6dd;padding:12px 0}svg{max-width:100%;height:auto;background:#fff}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid #ccd6dd;padding:6px;text-align:left}th{background:#e9f3f8}.note{background:#fff6d9;padding:12px}.back{color:#176b87;font-weight:700;text-decoration:none}</style></head><body><p><a class='back' href='../../../../index.html'>&larr; Back to reports home</a></p><h1>Low-resource staged parameter estimation</h1><p>Generated $(Dates.now()). Fits use <code>GrowthParameterEstimation</code> on the workbook-derived mean-cell trajectories. Points are observed values; lines are the BIC-selected ODE fit.</p><div class='note'><b>Identifiability note.</b> Low-resource workbooks do not provide a time-resolved resource concentration, and Ce1 has one treatment level. Consequently, the resource/treatment effect is a constant lineage-specific loss term; a Hill IC50 curve is not estimated. “Joint” Run 1 + Run 2 fits share shape parameters after each run is normalised at its observed day-zero value.</div><h2>Winning model ledger</h2><table><tr><th>Stage</th><th>Condition</th><th>Fit scope</th><th>Winning model</th><th>BIC</th><th>Inherited parameters</th></tr>$summary_rows</table><h2>Fits and BIC model comparisons</h2>$(join(panels,"\n"))</body></html>"""
open(joinpath(OUT,"report.html"),"w") do io; write(io,html); end
println("Wrote $(joinpath(OUT,"report.html")) and bic_model_ranking.csv")
