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
function delayed_death_logistic!(du, u, p, t)
    r, K, d, tdeath = p
    du[1] = r*u[1]*(1-u[1]/max(K, 1e-8)) - (t >= tdeath ? d : 0.0)*u[1]
end
function smooth_delayed_death_logistic!(du, u, p, t)
    r, K, d, tdeath, width = p
    onset = 1 / (1 + exp(-clamp((t - tdeath) / max(width, 1e-4), -30.0, 30.0)))
    du[1] = r*u[1]*(1-u[1]/max(K, 1e-8)) - d*onset*u[1]
end
function transit_delayed_kill!(du, u, p, t)
    r, K, kdamage, ktransit, kdeath = p
    N, D1, D2 = u
    growth = r*N*(1-N/max(K, 1e-8))
    du[1] = growth - kdamage*N
    du[2] = kdamage*N - ktransit*D1
    du[3] = ktransit*D1 - (ktransit + kdeath)*D2
end

function candidate_specs(scale)
    cap = max(2.0, scale * 1.25)
    return Dict(
        "logistic" => (model=logistic!, p0=[0.35, cap], bounds=[(1e-6, 4.0), (1e-3, max(10.0, scale*20))], initial=y0->[y0]),
        "logistic_plus_loss" => (model=loss_logistic!, p0=[0.5, cap, 0.08], bounds=[(1e-6, 4.0), (1e-3, max(10.0, scale*20)), (0.0, 3.0)], initial=y0->[y0]),
        "logistic_with_lag" => (model=delay_logistic!, p0=[0.5, cap, 0.5], bounds=[(1e-6, 4.0), (1e-3, max(10.0, scale*20)), (0.0, 8.0)], initial=y0->[y0]),
        "logistic_with_delayed_death" => (model=delayed_death_logistic!, p0=[0.6, cap, 0.35, 7.0], bounds=[(1e-6,4.0),(1e-3,max(10.0,scale*20)),(0.0,3.0),(0.0,13.5)], initial=y0->[y0]),
        "logistic_with_smooth_delayed_death" => (model=smooth_delayed_death_logistic!, p0=[0.6, cap, 0.35, 7.0, 0.8], bounds=[(1e-6,4.0),(1e-3,max(10.0,scale*20)),(0.0,3.0),(0.0,13.5),(0.05,5.0)], initial=y0->[y0]),
        "transit_compartment_delayed_kill" => (model=transit_delayed_kill!, p0=[0.6,cap,0.1,0.4,0.25], bounds=[(1e-6,4.0),(1e-3,max(10.0,scale*20)),(0.0,3.0),(1e-4,5.0),(0.0,5.0)], initial=y0->[y0,0.0,0.0]),
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
            result = GPE.run_joint_fit(spec.model, [(x=x, y=y, state_index=1)], spec.initial(y[1]), spec.p0;
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
        result = GPE.run_joint_fit(spec.model, data, spec.initial(1.0), spec.p0; bounds=spec.bounds, maxiters=2_000, optimizer=:nelder_mead, reltol=1e-7, abstol=1e-7)
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
    height = max(height, 32 + 42 * length(rows))
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

all_rows = NamedTuple[]
run1_panels = String[]; run2_panels = String[]; joint_panels = String[]; coculture_panels = String[]
ratio_buckets = Dict{String,Vector{NamedTuple}}()

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
    ratio_match = match(r" ([13]-[13])$", g.sheet[1])
    if ratio_match !== nothing
        ratio = ratio_match.captures[1]
        density_match = match(r" (\d+k) [13]-[13]$", g.sheet[1])
        density = density_match === nothing ? "" : density_match.captures[1]
        lineage = replace(g.sheet[1], r"^[12] - " => "")
        lineage = replace(lineage, r" \d+k [13]-[13]$" => "")
        ratio_key = "$(run) | $(g.workbook[1]) | $(density) | $(ratio)"
        push!(get!(ratio_buckets, ratio_key, NamedTuple[]), (lineage=lineage, density=density, ratio=ratio, run=run, x=x, y=y, pred=pred, ranked=ranked))
    else
        panel = "<section><h3>$(esc(name)) — $(esc(run))</h3>$(svg_panel("Observed points and winning fit: $(best.model)", [(name=run,x=x,y=y)], [pred]))$(bic_svg(ranked))</section>"
        push!(run == "Run 2" ? run2_panels : run1_panels, panel)
    end
end

# Sheets sharing a seeding-ratio label (for example 1-1, 1-3, or 3-1) are
# paired sensitive/resistant measurements.  Keep their independently selected
# fits, but render observed values and both fitted curves together.
for entries in values(ratio_buckets)
    first_entry = entries[1]
    title = "$(first_entry.run) — $(first_entry.density) total — $(first_entry.ratio) sensitive:resistant ratio"
    series = [(name=e.lineage, x=e.x, y=e.y) for e in entries]
    predictions = [e.pred for e in entries]
    bic_charts = join(["<h4>$(esc(e.lineage)) model comparison</h4>$(bic_svg(e.ranked))" for e in entries], "")
    panel = "<section><h3>$(esc(title))</h3>$(svg_panel("Observed and fitted growth for both lineages", series, predictions))$(bic_charts)</section>"
    push!(first_entry.run == "Run 2" ? run2_panels : run1_panels, panel)
end

for (key, runs) in run_buckets
    length(runs) == 2 || continue
    ranked, fits = fit_joint(runs); best = ranked[1]
    pred = [[fits[best.model].predictions[i][j] * runs[i].y[1] for j in eachindex(runs[i].x)] for i in eachindex(runs)]
    append!(all_rows, [(stage="LR-1 monoculture/lineage", condition=key, fit_scope="Run 1 + Run 2 (joint normalized shape)", model=r.model, bic=r.bic, sse=r.sse, params=r.params, inherited="shared r/K/loss across runs; run-specific observed N0") for r in ranked])
    push!(joint_panels, "<section><h3>$(esc(key)) — joint Run 1 + Run 2</h3>$(svg_panel("Joint fitted shape restored to each run's observed N0", [(name=r.name,x=r.x,y=r.y) for r in runs], pred))$(bic_svg(ranked))</section>")
end

# Stage LR-2/LR-3: Ce0 provides untreated coculture baselines. Ce1 inherits those
# parameters and compares constant loss with delayed Hill killing at C = IC50 = 1 µM.
# With one concentration, Hill is fixed to 1; Emax and delayed-onset times remain fit.
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
    pred=fits[best.model].predictions; push!(coculture_panels,"<section><h3>$(esc(g.workbook[1])) — untreated coculture</h3>$(svg_panel("Observed and winning inherited-baseline candidate",[(name="Sensitive",x=x,y=S),(name="Resistant",x=x,y=R)],pred))$(bic_svg(rows))</section>")
end

for g in groupby(filter(row -> startswith(lowercase(row.workbook), "ce1"), ce), :workbook)
    ratio=replace(g.workbook[1],"ce1_"=>""); haskey(ce0,ratio) || continue
    x,S,R=ce_pair(g); base=ce0[ratio].fit.params; data=[(x=x,y=S,state_index=1),(x=x,y=R,state_index=2)]
    inherited_constant_fun = if ce0[ratio].model == "competition_logistic"
        (du,u,p,t)->begin rS,KS,aSR,rR,KR,aRS=base; dS,dR=p; du[1]=rS*u[1]*(1-(u[1]+aSR*u[2])/max(KS,1e-8))-dS*u[1]; du[2]=rR*u[2]*(1-(u[2]+aRS*u[1])/max(KR,1e-8))-dR*u[2] end
    else
        (du,u,p,t)->begin rS,KS,rR,KR=base; dS,dR=p; du[1]=rS*u[1]*(1-u[1]/max(KS,1e-8))-dS*u[1]; du[2]=rR*u[2]*(1-u[2]/max(KR,1e-8))-dR*u[2] end
    end
    inherited_delayed_hill_fun = if ce0[ratio].model == "competition_logistic"
        (du,u,p,t)->begin
            rS,KS,aSR,rR,KR,aRS=base; emaxS,emaxR,tlagS,tlagR=p
            effectS=emaxS*0.5*(t >= tlagS ? 1.0 : 0.0); effectR=emaxR*0.5*(t >= tlagR ? 1.0 : 0.0)
            du[1]=rS*u[1]*(1-(u[1]+aSR*u[2])/max(KS,1e-8))-effectS*u[1]
            du[2]=rR*u[2]*(1-(u[2]+aRS*u[1])/max(KR,1e-8))-effectR*u[2]
        end
    else
        (du,u,p,t)->begin
            rS,KS,rR,KR=base; emaxS,emaxR,tlagS,tlagR=p
            effectS=emaxS*0.5*(t >= tlagS ? 1.0 : 0.0); effectR=emaxR*0.5*(t >= tlagR ? 1.0 : 0.0)
            du[1]=rS*u[1]*(1-u[1]/max(KS,1e-8))-effectS*u[1]
            du[2]=rR*u[2]*(1-u[2]/max(KR,1e-8))-effectR*u[2]
        end
    end
    specs=Dict(
        "inherited_$(ce0[ratio].model)_plus_constant_treatment_loss" => (inherited_constant_fun,[0.1,0.1],[(0.0,3.0),(0.0,3.0)]),
        "inherited_$(ce0[ratio].model)_plus_delayed_Hill_kill_IC50_1uM" => (inherited_delayed_hill_fun,[0.3,0.3,7.0,7.0],[(0.0,6.0),(0.0,6.0),(0.0,13.5),(0.0,13.5)]),
    )
    rows=NamedTuple[]; fits=Dict{String,Any}()
    for (model,(fun,p0,bounds)) in specs
        f=GPE.run_joint_fit(fun,data,[S[1],R[1]],p0;bounds=bounds,maxiters=3000,optimizer=:nelder_mead,reltol=1e-7,abstol=1e-7)
        fits[model]=f; push!(rows,(model=model,bic=f.bic,sse=f.sse,params=join(round.(f.params;sigdigits=5),";")))
    end
    sort!(rows; by = row -> row.bic); best=rows[1]
    append!(all_rows,[(stage="LR-3 treated coculture",condition=g.workbook[1],fit_scope="available ratio",model=r.model,bic=r.bic,sse=r.sse,params=r.params,inherited="$(ce0[ratio].model) parameters from ce0_$(ratio); C=IC50=1 µM, Hill=1; fitted treatment terms") for r in rows])
    push!(coculture_panels,"<section><h3>$(esc(g.workbook[1])) — treated coculture</h3>$(svg_panel("Observed and winning inherited fit (C = IC50 = 1 µM)",[(name="Sensitive",x=x,y=S),(name="Resistant",x=x,y=R)],fits[best.model].predictions))$(bic_svg(rows))</section>")
end

results=DataFrame(all_rows); CSV.write(joinpath(OUT,"bic_model_ranking.csv"),results)
summary_rows = join(["<tr><td>$(esc(r.stage))</td><td>$(esc(r.condition))</td><td>$(esc(r.fit_scope))</td><td>$(esc(r.model))</td><td>$(round(r.bic;digits=2))</td><td>$(esc(r.inherited))</td></tr>" for r in eachrow(results) if r.bic == minimum(results.bic[(results.stage .== r.stage) .& (results.condition .== r.condition) .& (results.fit_scope .== r.fit_scope)])], "\n")
tab_panels = """<div class='tabs' role='tablist' aria-label='Fit comparison groups'><button class='tab active' role='tab' aria-selected='true' data-tab='run1'>Run 1</button><button class='tab' role='tab' aria-selected='false' data-tab='run2'>Run 2</button><button class='tab' role='tab' aria-selected='false' data-tab='joint'>Joint</button><button class='tab' role='tab' aria-selected='false' data-tab='coculture'>Coculture</button></div><div class='tab-panel active' id='run1' role='tabpanel'>$(join(run1_panels,"\n"))</div><div class='tab-panel' id='run2' role='tabpanel'>$(join(run2_panels,"\n"))</div><div class='tab-panel' id='joint' role='tabpanel'>$(join(joint_panels,"\n"))</div><div class='tab-panel' id='coculture' role='tabpanel'>$(join(coculture_panels,"\n"))</div>"""
html = """<!doctype html><html><head><meta charset='utf-8'><title>Low-resource staged parameter estimation</title><style>body{font:14px system-ui;margin:32px;color:#17212b;max-width:1250px}h1{color:#124b6e}section{border-top:1px solid #ccd6dd;padding:12px 0}svg{max-width:100%;height:auto;background:#fff}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid #ccd6dd;padding:6px;text-align:left}th{background:#e9f3f8}.note{background:#fff6d9;padding:12px}.back{color:#176b87;font-weight:700;text-decoration:none}.tabs{display:flex;gap:6px;flex-wrap:wrap;margin:12px 0}.tab{border:1px solid #176b87;background:#f7fbfd;color:#124b6e;border-radius:5px;padding:8px 14px;font-weight:700;cursor:pointer}.tab.active,.tab:hover{background:#176b87;color:#fff}.tab-panel{display:none}.tab-panel.active{display:block}</style></head><body><p><a class='back' href='../../../../index.html'>&larr; Back to reports home</a></p><h1>Low-resource staged parameter estimation</h1><p>Generated $(Dates.now()). Fits use <code>GrowthParameterEstimation</code> on the workbook-derived mean-cell trajectories. Points are observed values; lines are the BIC-selected ODE fit.</p><div class='note'><b>Late-drop and treatment assumptions.</b> The model set retains the original logistic, constant-loss, and delayed-growth candidates and adds abrupt delayed death, smooth delayed death, and a transit-compartment delayed-kill model. For Ce1, drug concentration and IC50 are both fixed at 1 µM, so the Hill effect is 0.5 × Emax with Hill fixed at 1; constant and delayed-kill alternatives are compared. “Joint” Run 1 + Run 2 fits share shape parameters after each run is normalised at its observed day-zero value.</div><h2>Winning model ledger</h2><table><tr><th>Stage</th><th>Condition</th><th>Fit scope</th><th>Winning model</th><th>BIC</th><th>Inherited parameters</th></tr>$summary_rows</table><h2>Fits and BIC model comparisons</h2>$tab_panels<script>for(const b of document.querySelectorAll('.tab'))b.addEventListener('click',()=>{document.querySelectorAll('.tab,.tab-panel').forEach(x=>{x.classList.remove('active');if(x.classList.contains('tab'))x.setAttribute('aria-selected','false')});b.classList.add('active');b.setAttribute('aria-selected','true');document.getElementById(b.dataset.tab).classList.add('active')})</script></body></html>"""
open(joinpath(OUT,"report.html"),"w") do io; write(io,html); end
println("Wrote $(joinpath(OUT,"report.html")) and bic_model_ranking.csv")
