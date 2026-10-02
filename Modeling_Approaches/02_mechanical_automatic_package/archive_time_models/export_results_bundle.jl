# Writes docs/outputs/results/bundle.js (window.RESULTS = {...}): everything the Summary,
# Presentation and Explorer pages draw — data, the bootstrap-optimal model of every
# fit group, their equations, parameters, bootstrap intervals, curves and bands, plus what the
# in-browser simulator needs to re-run each model. Included at the end of run_bootstrap_fits.jl,
# so every full run refreshes the site.

const RESULTS_DIR = joinpath(ROOT, "docs", "outputs", "results")

# ── minimal JSON writer ──────────────────────────────────────────────────────
jnum(x) = isfinite(x) ? @sprintf("%.6g", x) : "null"
function jsonv(x)
    x === nothing && return "null"
    x isa Bool && return x ? "true" : "false"
    x isa Real && return jnum(Float64(x))
    x isa Symbol && return jsonv(string(x))
    if x isa AbstractString
        io = IOBuffer(); print(io, '"')
        for c in x
            c == '"' ? print(io, "\\\"") : c == '\\' ? print(io, "\\\\") : c == '\n' ? print(io, "\\n") :
            c < ' ' ? print(io, "\\u", string(UInt16(c); base=16, pad=4)) : print(io, c)
        end
        print(io, '"'); return replace(String(take!(io)), "</" => "<\\/")
    end
    x isa NamedTuple && return jsonv(Dict(string(k) => v for (k, v) in pairs(x)))
    x isa AbstractDict && return "{" * join(["$(jsonv(string(k))):$(jsonv(v))" for (k, v) in x], ",") * "}"
    x isa Union{AbstractVector,Tuple} && return "[" * join(jsonv.(collect(x)), ",") * "]"
    x isa AbstractMatrix && return jsonv([x[:, j] for j in axes(x, 2)])
    error("jsonv: unsupported $(typeof(x))")
end

# ── symbolic equations (no numbers) ──────────────────────────────────────────
function sym_tex(g, c)
    if g.kind == :s1
        X = stage1[g.line].sym; G = "r_{$X}\\,$X\\left(1-\\frac{$X}{K_{$X}}\\right)"
        c.name == "logistic" && return ["\\frac{d$X}{dt}=$G"]
        c.name == "theta_logistic" && return ["\\frac{d$X}{dt}=r_{$X}\\,$X\\left[1-\\left(\\frac{$X}{K_{$X}}\\right)^{\\theta_{$X}}\\right]"]
        return ["\\frac{d$X}{dt}=r_{$X}\\,$X\\,\\log\\left(\\frac{K_{$X}}{$X}\\right)"]
    elseif g.kind == :s1b
        return [c.tex]
    elseif g.kind == :s2
        X = stage1[g.line].sym; G = "r_{$X}\\,$X\\left(1-\\frac{$X}{K_{$X}}\\right)"; A = "\\left(1-e^{-\\lambda_{$X} t}\\right)"
        H = "\\frac{D^{h_{$X}}}{EC_{50}^{h_{$X}}+D^{h_{$X}}}"
        c.name == "delayed_hill_kill" && return ["\\frac{d$X}{dt}=$G-\\beta_{$X}$A$H\\,$X"]
        c.name == "delayed_linear_kill" && return ["\\frac{d$X}{dt}=$G-\\beta_{$X}$A\\,D\\,$X"]
        c.name == "immediate_hill_kill" && return ["\\frac{d$X}{dt}=$G-\\beta_{$X}$H\\,$X"]
        return ["\\frac{d$X}{dt}=$G-\\beta_{$X}\\,D\\,$X"]
    elseif g.kind == :s3
        inh = c.name in S3_INHERITED_KILL
        eq3(X, Y, a) = begin
            crowd = "α_$a" in c.pnames ? "$X+\\alpha_{$a}$Y" : "$X+$Y"
            s = "\\frac{d$X}{dt}=r_{$X}$X\\left(1-\\frac{$crowd}{K}\\right)"
            "d_$X" in c.pnames && (s *= "-d_{$X}$X")
            (inh || "β_$X" in c.pnames) && (s *= "-\\mathbf{1}_{\\mathrm{Ce1}}\\beta_{$X}A_{$X}(t)H_{$X}(D)\\,$X")
            s
        end
        return [eq3("N", "C", "NC"), eq3("C", "N", "CN")]
    else
        d = s3b_defaults(c)
        eq3b(X, Y, a, dv, av) = begin
            crowd = av == 0.0 && !("α_$a" in c.pnames) ? "$X" : av == 1.0 && !("α_$a" in c.pnames) ? "$X+$Y" : "$X+\\alpha_{$a}$Y"
            s = "\\frac{d$X}{dt}=r_{$X}$X\\left(1-\\frac{$crowd}{K}\\right)"
            (dv != 0.0 || "d_$(X)^LR" in c.pnames) && (s *= "-d^{LR}_{$X}u_R\\,$X")
            s
        end
        return [eq3b("N", "C", "NC", d["dN"], d["aNC"]), eq3b("C", "N", "CN", d["dC"], d["aCN"])]
    end
end

# Stage 3b: values used when a term is not fitted
function s3b_defaults(c)
    n = c.name
    n == "lr_cc_sharedK"       && return Dict("K" => K3, "aNC" => 1.0, "aCN" => 1.0, "dN" => dN1b, "dC" => dC1b)
    n == "lr_cc_crowd_no_loss" && return Dict("K" => K3, "aNC" => aNC3, "aCN" => aCN3, "dN" => 0.0, "dC" => 0.0)
    n == "lr_cc_alphaCN_only"  && return Dict("K" => K3, "aNC" => 0.0, "aCN" => aCN3, "dN" => 0.0, "dC" => 0.0)
    Dict("K" => K3, "aNC" => aNC3, "aCN" => aCN3, "dN" => dN1b, "dC" => dC1b)
end

# inherited / fixed values shown next to the fitted ones
function inherited_rows(g, c)
    if g.kind == :s1
        return []
    elseif g.kind == :s1b
        s1 = stage1[g.line]; X = s1.sym
        return [(sym="r_{$X}", value=s1.r, note="Stage 1"), (sym="K_{$X}", value=s1.K, note="Stage 1"), (sym="u_R", value=U_R, note="0.5% FBS")]
    elseif g.kind == :s2
        s1 = stage1[g.line]; X = s1.sym
        rows = [(sym="r_{$X}", value=s1.r, note="Stage 1"), (sym="K_{$X}", value=s1.K, note="Stage 1")]
        occursin("hill", c.name) && push!(rows, (sym="EC_{50}", value=EC50, note="fixed, µM"))
        return rows
    elseif g.kind == :s3
        rows = [(sym="r_N", value=rN, note="Stage 1"), (sym="r_C", value=rC, note="Stage 1"), (sym="\\lambda_N", value=λN, note="Stage 2"),
                (sym="\\lambda_C", value=λC, note="Stage 2"), (sym="H_N(1\\,\\mu M)", value=HN, note="Stage 2"), (sym="H_C(1\\,\\mu M)", value=HC, note="Stage 2")]
        c.name in S3_INHERITED_KILL && append!(rows, [(sym="\\beta_N", value=βN2, note="Stage 2"), (sym="\\beta_C", value=βC2, note="Stage 2")])
        "α_NC" in c.pnames || c.name in S3_INHERITED_KILL || push!(rows, (sym="\\alpha_{NC}=\\alpha_{CN}", value=1.0, note="one-for-one"))
        return rows
    else
        d = s3b_defaults(c); rows = [(sym="r_N", value=rN, note="Stage 1"), (sym="r_C", value=rC, note="Stage 1"), (sym="u_R", value=U_R, note="0.5% FBS")]
        for (k, sym, src) in (("K", "K", "Stage 3"), ("aNC", "\\alpha_{NC}", "Stage 3"), ("aCN", "\\alpha_{CN}", "Stage 3"),
                              ("dN", "d^{LR}_N", "Stage 1b"), ("dC", "d^{LR}_C", "Stage 1b"))
            name = Dict("K" => "K", "aNC" => "α_NC", "aCN" => "α_CN", "dN" => "d_N^LR", "dC" => "d_C^LR")[k]
            name in c.pnames && continue
            push!(rows, (sym=sym, value=d[k], note=d[k] == 0.0 ? "term removed" : (k in ("aNC", "aCN") && d[k] == 1.0) ? "one-for-one" : src))
        end
        return rows
    end
end

# TeX for a fitted parameter name as used in the Julia code
function ptex(n)
    m = Dict("θ" => "\\theta", "ρ" => "\\rho", "β" => "\\beta", "λ" => "\\lambda", "α" => "\\alpha", "δ" => "\\delta", "γ" => "\\gamma", "κ" => "\\kappa", "τ" => "\\tau")
    base, rest = n, ""
    if occursin("^", n)
        base, sup = split(n, "^"; limit=2); rest = "^{$sup}"
    end
    parts = split(base, "_"; limit=2)
    b = get(m, string(parts[1]), string(parts[1]))
    length(parts) == 2 ? "$(b)_{$(parts[2])}$rest" : "$b$rest"
end

# what the in-browser simulator needs to re-run a model
function js_spec(g, c)
    if g.kind == :s1
        return (id=c.name, consts=Dict{String,Any}(), u0=cu0(c, g))
    elseif g.kind == :s1b
        s1 = stage1[g.line]; return (id=c.name, consts=Dict{String,Any}("r" => s1.r, "K" => s1.K, "uR" => U_R), u0=cu0(c, g))
    elseif g.kind == :s2
        s1 = stage1[g.line]; return (id=c.name, consts=Dict{String,Any}("r" => s1.r, "K" => s1.K, "EC50" => EC50, "doses" => stage2[g.line].doses), u0=cu0(c, g))
    elseif g.kind == :s3
        inh = c.name in S3_INHERITED_KILL
        return (id="coculture", consts=Dict{String,Any}("rN" => rN, "rC" => rC, "lN" => λN, "lC" => λC, "HN" => HN, "HC" => HC,
            "betaN" => inh ? βN2 : 0.0, "betaC" => inh ? βC2 : 0.0, "envDose" => ENV_DOSE, "pnames" => c.pnames), u0=cu0(c, g))
    else
        return (id="lr_coculture", consts=Dict{String,Any}("rN" => rN, "rC" => rC, "uR" => U_R, "defaults" => s3b_defaults(c), "pnames" => c.pnames), u0=cu0(c, g))
    end
end

# data series: what is plotted, and which model state each series is compared with
function data_series(g)
    if g.kind in (:s1, :s1b)
        s = g.kind == :s1 ? stage1[g.line] : stage1b[g.line]
        bad = g.kind == :s1b ? s.bad : falses(length(s.x))
        return [(panel=1, name="observed", color=BLUE, filled=true, dashed=false, state=1, x=s.x[.!bad], y=s.y[.!bad], xbad=s.x[bad], ybad=s.y[bad])]
    elseif g.kind == :s2
        cols = [PURPLE, GREEN, ORANGE]
        return [(panel=1, name="$(dose) µM", color=cols[i], filled=true, dashed=false, state=i, x=d.x, y=d.y, xbad=Float64[], ybad=Float64[])
                for (i, (dose, d)) in enumerate(zip(stage2[g.line].doses, g.datasets))]
    elseif g.kind == :s3
        out = NamedTuple[]
        for (pn, ratio) in enumerate(["25S75R", "50S50R", "75S25R"])
            i0 = findfirst(e -> e.ratio == ratio && e.cond == "Ce0", coculture_layout); i1 = findfirst(e -> e.ratio == ratio && e.cond == "Ce1", coculture_layout)
            d0, d1 = co_data[i0], co_data[i1]
            append!(out, [(panel=pn, name="Ce0 sensitive", color=BLUE, filled=true, dashed=false, state=2i0-1, x=d0.S[1], y=d0.S[2], xbad=Float64[], ybad=Float64[]),
                          (panel=pn, name="Ce0 resistant", color=ORANGE, filled=true, dashed=false, state=2i0, x=d0.R[1], y=d0.R[2], xbad=Float64[], ybad=Float64[]),
                          (panel=pn, name="Ce1 sensitive", color=BLUE, filled=false, dashed=true, state=2i1-1, x=d1.S[1], y=d1.S[2], xbad=Float64[], ybad=Float64[]),
                          (panel=pn, name="Ce1 resistant", color=ORANGE, filled=false, dashed=true, state=2i1, x=d1.R[1], y=d1.R[2], xbad=Float64[], ybad=Float64[])])
        end
        return out
    else
        out = NamedTuple[]
        for (i, c) in enumerate(lr_cc)
            push!(out, (panel=i, name="sensitive", color=BLUE, filled=true, dashed=false, state=2i-1, x=c.xN[.!c.bN], y=c.yN[.!c.bN], xbad=c.xN[c.bN], ybad=c.yN[c.bN]))
            push!(out, (panel=i, name="resistant", color=ORANGE, filled=true, dashed=false, state=2i, x=c.xC[.!c.bC], y=c.yC[.!c.bC], xbad=c.xC[c.bC], ybad=c.yC[c.bC]))
        end
        return out
    end
end

panel_titles(g) = g.kind == :s3 ? ["25:75 sensitive:resistant", "50:50", "75:25"] :
                  g.kind == :s3b ? ["Run $(c.run), $(c.ratio) S:R" for c in lr_cc] : [""]
main_panel(g) = g.kind == :s3 ? 2 : g.kind == :s3b ? 2 : 1

const GROUP_TEXT = Dict(
    :s1 => "Grown alone at 30k in 10% FBS. r and K from this stage are carried into every later stage.",
    :s1b => "Grown alone at 0.5% FBS (30k sheet). r and K are inherited from Stage 1; only the low-resource terms are fitted.",
    :s2 => "Cisplatin at 0.62, 1.0 and 1.47 µM, fitted jointly. r and K inherited from Stage 1; EC50 fixed at 1 µM.",
    :s3 => "Untreated (Ce0) and 1 µM treated (Ce1) cocultures at three ratios, all fitted at once with one shared carrying capacity K.",
    :s3b => "Low-resource cocultures (0.5% FBS, runs 1 and 2, three ratios), all fitted at once with one shared carrying capacity K.")

function model_entry(g, r, j, role; band=nothing)
    c = r.cands[j]; tmax = g.kind == :s1b ? 20.0 : XMAX
    ts, sol = dense(c.model, cu0(c, g), c.params, tmax; dt=0.1)
    states = sort(unique(s.state for s in data_series(g)))
    params = [begin
                  l, m, h = qs(r.pars[j][:, i])
                  (name=n, tex=ptex(n), value=c.params[i], lo=l, hi=h, removable=removable(l, h, c.bounds[i]),
                   atbound=isapprox(c.params[i], c.bounds[i][2]; rtol=1e-3))
              end for (i, n) in enumerate(c.pnames)]
    e = Dict{String,Any}("role" => role, "name" => c.name, "label" => c.label, "k" => r.ks[j], "bic" => c.bic,
        "dbic" => c.bic - minimum(x.bic for x in r.cands), "pick" => pct(r.pick, j), "lowest" => pct(r.winner, j),
        "tex" => sym_tex(g, c), "params" => params, "inherited" => inherited_rows(g, c),
        "curve" => Dict("t" => ts, "y" => Dict(string(s) => [u[s] for u in sol.u] for s in states)),
        "sim" => js_spec(g, c), "bounds" => [collect(b) for b in c.bounds])
    if g.kind == :s1b
        e["day30"] = c.x30; e["reachesZero"] = c.extinct
    end
    if band !== nothing
        e["band"] = Dict("t" => band.ts, "lo" => Dict(string(s) => band.lo[:, s] for s in states), "hi" => Dict(string(s) => band.hi[:, s] for s in states))
    end
    e
end

function export_results_bundle(groups, results; B=BOOT_B)
    out = Any[]
    for g in groups
        r = results[g.id].res
        jo = r.selected
        gd = Dict{String,Any}("id" => g.id, "stage" => g.stage, "title" => g.title, "kind" => string(g.kind), "text" => GROUP_TEXT[g.kind],
            "series" => data_series(g), "panels" => panel_titles(g), "mainPanel" => main_panel(g),
            "optimal" => model_entry(g, r, jo, "optimal"; band=results[g.id].band), "candidates" => length(r.cands),
            "bicBest" => Dict("label" => r.cands[1].label, "k" => r.ks[1], "bic" => r.cands[1].bic, "isOptimal" => jo == 1),
            "xmax" => g.kind == :s1b ? 20.0 : XMAX,
            "ymax" => g.kind == :s1b ? ceil(maximum(maximum(s.y) for s in data_series(g)) * 1.2; sigdigits=2) : YMAX)
        push!(out, gd)
    end
    bundle = Dict{String,Any}("generated" => Dates.format(now(), "yyyy-mm-dd HH:MM"), "replicates" => B,
        "perParamBIC" => PER_PARAM_BIC, "groups" => out)
    mkpath(RESULTS_DIR)
    write(joinpath(RESULTS_DIR, "bundle.js"), "window.RESULTS = " * jsonv(bundle) * ";\n")
    println("Wrote $(joinpath(RESULTS_DIR, "bundle.js"))")
end
