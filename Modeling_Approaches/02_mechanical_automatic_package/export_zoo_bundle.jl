# Writes docs/outputs/results/bundle.js (window.RESULTS) from run_zoo_fits.jl: data, the
# optimal model (fewest free parameters within 2 BIC of the lowest) of every fit group with its equations, parameters, bootstrap
# intervals, curves and band, the full candidate ranking, and the system spec the in-browser
# simulator (site.js, MODELS.zoo) needs to re-run it.

const RESULTS_DIR = joinpath(ROOT, "docs", "outputs", "results")

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

line_syms(g) = g.kind in (:s3, :s3b, :s4) ? ["N", "C"] : [g.line == "A2780" ? "N" : "C"]

function crowd_tex(a::Slot, X, Y, name)
    (a.i == 0 && a.v == 0) && return X
    (a.i == 0 && a.v == 1) && return "$X+$Y"
    "$X+\\alpha_{$name}$Y"
end

function zoo_tex(g, c)
    S = c.sys
    if g.kind in (:s3, :s3b, :s4)
        return vcat(line_tex(S.lines[1], "N", crowd_tex(S.alpha[1][2], "N", "C", "NC"), "K"),
                    line_tex(S.lines[2], "C", crowd_tex(S.alpha[2][1], "C", "N", "CN"), "K"))
    end
    X = line_syms(g)[1]
    line_tex(S.lines[1], X, X, "K_{$X}")
end

# values that are carried over or fixed, shown next to the fitted ones
function zoo_inherited(g, c)
    g.kind == :s1 && return NamedTuple[]
    g.kind == :s4 && return s4_inherited(c)
    rows = NamedTuple[]
    for (l, X) in enumerate(line_syms(g))
        L = c.sys.lines[l]
        push!(rows, (sym="r_{$X}", value=L.s.r.v, note="Stage 1"))
        g.kind in (:s1b, :s2) && push!(rows, (sym="K_{$X}", value=L.s.K.v, note="Stage 1"))
        L.growth == :theta && push!(rows, (sym="\\theta_{$X}", value=L.s.θ.v, note="Stage 1"))
        if g.kind in (:s3, :s3b)
            mech = c.mech[l]; src = g.kind == :s3 ? "Stage 2" : "Stage 1b"
            for k in mech.slots
                getfield(L.s, k).i > 0 && continue
                push!(rows, (sym=ptex(slot_name(k, X, L.input; combo=(:d0 in mech.slots))), value=getfield(L.s, k).v, note=src))
            end
            c.noinput && push!(rows, (sym="\\varepsilon_{$X}", value=0.0, note="low-resource effect off"))
        end
    end
    if g.kind in (:s3, :s3b)
        c.sys.lines[1].s.K.i == 0 && push!(rows, (sym="K", value=c.sys.lines[1].s.K.v, note="Stage 3"))
        for (a, nm) in ((c.sys.alpha[1][2], "\\alpha_{NC}"), (c.sys.alpha[2][1], "\\alpha_{CN}"))
            a.i == 0 && push!(rows, (sym=nm, value=a.v, note=a.v == 1 ? "one-for-one" : a.v == 0 ? "term removed" : "Stage 3"))
        end
    end
    g.kind in (:s1b, :s3b) && push!(rows, (sym="u_R", value=U_R, note="0.5% FBS"))
    if g.kind in (:s2, :s3)
        for (L, X) in zip(c.sys.lines, line_syms(g))
            L.action == :hill && push!(rows, (sym="EC_{50,$X}", value=ec50_of(L.s, Float64[], Float64), note="measured IC50, µM"))
        end
    end
    g.kind == :s3 && push!(rows, (sym="D", value=COCULTURE_DOSE, note="µM in Ce1 (0 in Ce0)"))
    rows
end

function zoo_sim(c)
    S = c.sys
    modict(M) = Dict{String,Any}("growth" => string(M.growth), "kernel" => string(M.kernel), "action" => string(M.action),
                                 "input" => string(M.input), "s" => Dict(string(k) => [getfield(M.s, k).i, getfield(M.s, k).v] for k in SLOT_KEYS if hasproperty(M.s, k)))
    lines = [begin
                 d = modict(L)
                 hasproperty(L, :extra) && (d["extra"] = [modict(M) for M in L.extra])     # extra mechanism modules (Stage 4)
                 hasproperty(L, :chi) && (d["chi"] = [L.chi.i, L.chi.v])                   # drug × low-resource interaction
                 ns_on(L) && (d["ns"] = true)                                              # Norton–Simon: drug losses × growth factor
                 d
             end for L in S.lines]
    envs = [Dict{String,Any}("D" => e.D, "R" => e.R, "lines" => e.lines .- 1) for e in S.envs]
    alpha = [[[a.i, a.v] for a in row] for row in S.alpha]
    Dict{String,Any}("id" => "zoo", "sys" => Dict{String,Any}("lines" => lines, "envs" => envs, "alpha" => alpha, "xi" => S.xi, "hb" => S.hb, "n" => S.n),
                     "u0" => c.u0, "clamp" => sys_clamp(S), "uR" => U_R, "EC50" => EC50)
end

function zoo_series(g)
    c = g.cands[1]; obs = sys_obs(c.sys)
    if g.kind == :s4
        return [(panel=1, name="sensitive", color=BLUE, filled=true, dashed=false, state=obs[1], x=combo.xN[.!combo.bN], y=combo.yN[.!combo.bN], xbad=combo.xN[combo.bN], ybad=combo.yN[combo.bN]),
                (panel=1, name="resistant", color=ORANGE, filled=true, dashed=false, state=obs[2], x=combo.xC[.!combo.bC], y=combo.yC[.!combo.bC], xbad=combo.xC[combo.bC], ybad=combo.yC[combo.bC])]
    end
    if g.kind == :s1
        d = untreated[g.line]
        return [(panel=1, name="observed", color=BLUE, filled=true, dashed=false, state=1, x=d.x, y=d.y, xbad=Float64[], ybad=Float64[])]
    elseif g.kind == :s1b
        m = lrmono[g.line]
        return [(panel=1, name="observed", color=BLUE, filled=true, dashed=false, state=1, x=m.x[.!m.bad], y=m.y[.!m.bad], xbad=m.x[m.bad], ybad=m.y[m.bad])]
    elseif g.kind == :s2
        t = treated[g.line]; cols = [PURPLE, GREEN, ORANGE]
        return [(panel=1, name="$(dose) µM", color=cols[i], filled=true, dashed=false, state=obs[i], x=x, y=y, xbad=Float64[], ybad=Float64[])
                for (i, (dose, (x, y))) in enumerate(zip(t.doses, t.data))]
    elseif g.kind == :s3
        out = NamedTuple[]
        for (pn, ratio) in enumerate(["25S75R", "50S50R", "75S25R"])
            i0 = findfirst(e -> e.ratio == ratio && e.cond == "Ce0", COCULTURE_LAYOUT); i1 = findfirst(e -> e.ratio == ratio && e.cond == "Ce1", COCULTURE_LAYOUT)
            d0, d1 = co_data[i0], co_data[i1]
            append!(out, [(panel=pn, name="Ce0 sensitive", color=BLUE, filled=true, dashed=false, state=obs[2i0-1], x=d0.S[1], y=d0.S[2], xbad=Float64[], ybad=Float64[]),
                          (panel=pn, name="Ce0 resistant", color=ORANGE, filled=true, dashed=false, state=obs[2i0], x=d0.R[1], y=d0.R[2], xbad=Float64[], ybad=Float64[]),
                          (panel=pn, name="Ce1 sensitive", color=BLUE, filled=false, dashed=true, state=obs[2i1-1], x=d1.S[1], y=d1.S[2], xbad=Float64[], ybad=Float64[]),
                          (panel=pn, name="Ce1 resistant", color=ORANGE, filled=false, dashed=true, state=obs[2i1], x=d1.R[1], y=d1.R[2], xbad=Float64[], ybad=Float64[])])
        end
        return out
    else
        out = NamedTuple[]
        for (i, c) in enumerate(lr_cc)
            push!(out, (panel=i, name="sensitive", color=BLUE, filled=true, dashed=false, state=obs[2i-1], x=c.xN[.!c.bN], y=c.yN[.!c.bN], xbad=c.xN[c.bN], ybad=c.yN[c.bN]))
            push!(out, (panel=i, name="resistant", color=ORANGE, filled=true, dashed=false, state=obs[2i], x=c.xC[.!c.bC], y=c.yC[.!c.bC], xbad=c.xC[c.bC], ybad=c.yC[c.bC]))
        end
        return out
    end
end

zoo_panels(g) = g.kind == :s4 ? ["1:1, low resource + 1 µM cisplatin"] : g.kind == :s3 ? ["25:75 sensitive:resistant", "50:50", "75:25"] :
                g.kind == :s3b ? ["Run $(c.run), $(c.ratio) S:R" for c in lr_cc] : [""]
zoo_main_panel(g) = g.kind in (:s3, :s3b) ? 2 : 1

const ZOO_GROUP_TEXT = Dict(
    :s1 => "Grown alone at 30k in 10% FBS. The selected growth law (fewest terms within 2 BIC) and its r, K (and θ) are carried into every later stage.",
    :s1b => "Grown alone at 0.5% FBS. Growth is carried over from Stage 1; only the low-resource mechanism is fitted, and every mechanism acts through states driven by u_R, never through time.",
    :s2 => "Cisplatin at 0.62, 1.0 and 1.47 µM, fitted jointly. Growth is carried over from Stage 1; delays and accumulation come from drug-effect or damaged-cell states driven by D, never from time. A Hill kill uses each line's measured IC50 as EC50 (1 and 8.34 µM); with the Norton–Simon option the kill is scaled by the line's growth factor, so crowded cells die less.",
    :s3 => "Untreated (Ce0) and 1 µM treated (Ce1) cocultures at three ratios, all fitted at once with one shared carrying capacity K. Each line keeps its Stage 2 drug mechanism.",
    :s3b => "Low-resource cocultures (0.5% FBS, runs 1 and 2, three ratios), all fitted at once with one shared carrying capacity K. Each line keeps its Stage 1b low-resource mechanism.",
    :s4 => "The combination culture: low resource (0.5% FBS) plus 1 µM cisplatin, 1:1 at 30k. Each line carries both its Stage 3 drug model and its Stage 3b low-resource mechanism; everything is inherited and candidates add at most a drug × low-resource interaction χ, a refitted K or low-resource strengths. The starting counts are estimated, because day 1 is more than 3-fold above the nominal seeding of 50 cells per line.")

function zoo_entry(g, r, j, band; curve=true)
    c = g.cands[j]; tmax = g.kind == :s1b ? 20.0 : XMAX
    ts, sol = curve ? dense(c.model, c.u0, c.params, tmax; dt=0.1) : (Float64[], (u=Vector{Float64}[],))
    states = sort(unique(s.state for s in zoo_series(g)))
    hasb = size(r.pars[j], 1) > 0
    params = [begin
                  l, m, h = hasb ? qs(r.pars[j][:, i]) : (NaN, NaN, NaN)
                  (name=n, tex=ptex(n), value=c.params[i], lo=l, hi=h, removable=hasb && removable(l, h, c.bounds[i]),
                   atbound=at_bound(c.params[i], c.bounds[i]))
              end for (i, n) in enumerate(c.pnames)]
    lowest = hasb ? pct(r.winner, j) : NaN; chosen = hasb ? pct(r.pick, j) : NaN
    e = Dict{String,Any}("role" => j == 1 ? "optimal" : "candidate", "name" => c.name, "label" => c.label, "family" => c.family, "k" => length(c.params),
        "bic" => c.bic, "dbic" => c.bic - g.cands[1].bic, "pick" => chosen, "lowest" => lowest,
        "tex" => zoo_tex(g, c), "params" => params, "inherited" => zoo_inherited(g, c),
        "curve" => curve ? Dict("t" => ts, "y" => Dict(string(s) => [u[s] for u in sol.u] for s in states)) : nothing,
        "sim" => zoo_sim(c), "bounds" => [collect(b) for b in c.bounds])
    g.kind == :s1b && hasproperty(c, :x30) && (e["day30"] = c.x30; e["reachesZero"] = c.extinct)
    if band !== nothing
        e["band"] = Dict("t" => band.ts, "lo" => Dict(string(s) => band.lo[:, s] for s in states), "hi" => Dict(string(s) => band.hi[:, s] for s in states))
    end
    e
end

# Stage 4: every fixed value of both modules of both lines, with its source
function s4_inherited(c)
    rows = NamedTuple[]
    for (l, X) in enumerate(["N", "C"])
        L = c.sys.lines[l]; M = L.extra[1]
        push!(rows, (sym="r_{$X}", value=L.s.r.v, note="Stage 1"))
        L.growth == :theta && push!(rows, (sym="\\theta_{$X}", value=L.s.θ.v, note="Stage 1"))
        for k in vcat(kernel_slots(L.kernel), action_slots(L.action), [:d])
            sl = getfield(L.s, k); sl.i > 0 && continue
            k == :d && sl.v == 0 && continue
            note = (k == :β && l == 2) ? (sl.v == c.βC2 ? "Stage 2" : "Stage 3") : k in (:h, :k) ? "Stage 2" : "Stage 3"
            push!(rows, (sym=ptex(slot_name(k, X, :D)), value=sl.v, note=note))
        end
        combo_ = M.s.d0.i > 0 || M.s.d0.v != 0
        for k in vcat(kernel_slots(M.kernel), action_slots(M.action), [:ε])
            sl = getfield(M.s, k); sl.i > 0 && continue
            k == :ε && sl.v == 1 && continue
            push!(rows, (sym=ptex(slot_name(k, X, :R; combo=combo_)), value=sl.v, note="Stage 1b / 3b"))
        end
    end
    c.sys.lines[1].s.K.i == 0 && push!(rows, (sym="K", value=c.sys.lines[1].s.K.v, note=c.crowd == :S3 ? "Stage 3" : "Stage 3b"))
    for (a, nm) in ((c.sys.alpha[1][2], "\\alpha_{NC}"), (c.sys.alpha[2][1], "\\alpha_{CN}"))
        a.i == 0 && push!(rows, (sym=nm, value=a.v, note=a.v == 1 ? "one-for-one" : c.crowd == :S3 ? "Stage 3" : "Stage 3b"))
    end
    push!(rows, (sym="u_R", value=U_R, note="0.5% FBS"), (sym="D", value=COMBO_DOSE, note="µM cisplatin"),
          [(sym="EC_{50,$X}", value=ec50_of(L.s, Float64[], Float64), note="measured IC50, µM") for (L, X) in zip(c.sys.lines, ("N", "C"))]...)
    rows
end

# what the Explorer's model builder needs to assemble any kernel × effect on a monoculture stage
function zoo_builder(g)
    g.kind in (:s1b, :s2) || return nothing
    c = g.cands[1]; L = c.sys.lines[1]; X = line_syms(g)[1]
    kern = [Dict("id" => string(k), "label" => lab, "slots" => string.(kernel_slots(k))) for (k, _, lab) in KERNELS]
    acts = L.input == :D ? [:loss, :hill, :damage] : [:loss, :damage, :capacity, :transit, :fade, :rho, :sep_r, :resource, :toxin]
    alab = Dict(:loss => L.input == :D ? "Linear kill" : "Loss", :hill => "Hill kill", :damage => "Damaged cells die", :capacity => "Lowers K",
                :transit => "Two-stage damage", :fade => "Slows growth + immediate loss", :rho => "Growth-rate and K multipliers",
                :sep_r => "Separate growth rate", :resource => "Resource depletion", :toxin => "Toxic by-product")
    aslots(a) = a == :toxin ? [:γ, :κ] : action_slots(a)
    slots = union(vcat([kernel_slots(k) for (k, _, _) in KERNELS]...), vcat([aslots(a) for a in acts]...), [:d0])
    Dict{String,Any}("growth" => string(L.growth), "input" => string(L.input), "r" => L.s.r.v, "K" => L.s.K.v, "θ" => L.s.θ.v,
        "envs" => [Dict("D" => e.D, "R" => e.R) for e in c.sys.envs], "X0" => [c.u0[i] for i in sys_obs(c.sys)],
        "kernels" => kern, "actions" => [Dict("id" => string(a), "label" => alab[a], "slots" => string.(aslots(a))) for a in acts],
        "slots" => Dict(string(k) => Dict("name" => slot_name(k, X, L.input), "tex" => ptex(slot_name(k, X, L.input)),
                                          "lo" => SLOT_BOUNDS[k][1], "hi" => SLOT_BOUNDS[k][2], "start" => SLOT_STARTS[k][1]) for k in slots),
        "uR" => U_R, "EC50" => ec50_of(L.s, Float64[], Float64),      # this line's EC50 (the builder makes one line)
        "ns" => ns_on(L))                                             # built drug models get the stage's Norton–Simon coupling
end

function selection_note(g)
    (isdefined(Main, :JOINT_INFO) && g.kind in (:s1b, :s2)) || return nothing
    if g.kind == :s2
        p = JOINT_INFO["s2"].pick
        return "Both lines share one cisplatin mechanism structure with their own rates (same parent line, same DNA-damage mechanism; resistance changes rates). Chosen on the joint BIC of the two lines: $(p.label), joint BIC $(fmt(p.bic)), k = $(p.k)."
    end
    i = JOINT_INFO["s1b"]
    "Shared low-resource structure tested: each line's own model (A, summed BIC $(fmt(i.A.bic)), k = $(i.A.k)) vs the best shared structure (B, $(i.B.label), BIC $(fmt(i.B.bic)), k = $(i.B.k), ΔBIC $(@sprintf("%+.2f", i.dbic))). " *
        (i.useB ? "B is within $(Int(RULE_BIC)) BIC with no more parameters, so both lines share it." : "B is not within $(Int(RULE_BIC)) BIC with no more parameters, so each line keeps its own model.")
end

# Align each plotted series with the fitted dataset of the same state, so the Explorer's in-browser SSE and BIC
# use exactly the fitted points: points the fit does not use (QC exclusions) move to the excluded markers, day 0
# (the seeding, which is the initial condition) is not plotted as data, and weighted fits attach each point's
# residual scale w.
function align_series(g, series)
    weighted = any(d -> hasproperty(d, :w), g.datasets)
    map(series) do s
        j = findfirst(d -> d.state_index == s.state, g.datasets)
        j === nothing && (@warn "no fitted dataset for series $(s.name) of $(g.id)"; return s)
        D = g.datasets[j]; pos = Dict(zip(D.x, eachindex(D.x))); wv = collect(Float64, dweights(D))
        used = [haskey(pos, x) for x in s.x]; other = .!used .& (s.x .> 0)
        out = merge(s, (x=s.x[used], y=s.y[used], xbad=vcat(s.xbad, s.x[other]), ybad=vcat(s.ybad, s.y[other])))
        weighted ? merge(out, (w=[wv[pos[t]] for t in out.x],)) : out
    end
end

# a parameter on either end of its search range (within 0.1% of the range): its interval is not data-driven
at_bound(v, (lo, hi)) = isfinite(v) && (abs(v - lo) <= 1e-3 * (hi - lo) || abs(v - hi) <= 1e-3 * (hi - lo))

function export_zoo_bundle(groups, results; B=BOOT_B)
    out = Any[]
    for g in groups
        r = results[g.id].res
        ranking = [Dict{String,Any}("name" => c.name, "label" => c.label, "family" => c.family, "k" => length(c.params), "bic" => c.bic,
                                    "dbic" => c.bic - g.cands[1].bic, "lowest" => size(r.bic, 1) > 0 ? pct(r.winner, j) : NaN,
                                    "pick" => size(r.bic, 1) > 0 ? pct(r.pick, j) : NaN, "compartments" => ncompartments(c.sys.lines[1].kernel, c.sys.lines[1].action))
                   for (j, c) in enumerate(g.cands)]
        series = align_series(g, zoo_series(g))
        push!(out, Dict{String,Any}("id" => g.id, "stage" => stage_label(g.stage), "title" => g.title, "kind" => string(g.kind), "text" => ZOO_GROUP_TEXT[g.kind],
            "series" => series, "panels" => zoo_panels(g), "mainPanel" => zoo_main_panel(g),
            "optimal" => zoo_entry(g, r, 1, results[g.id].band), "candidates" => length(g.cands), "ranking" => ranking, "selectionNote" => selection_note(g),
            "models" => [zoo_entry(g, r, j, nothing; curve=false) for j in eachindex(g.cands)], "builder" => zoo_builder(g),
            "bicBest" => (jb = argmin([c.bic for c in g.cands]); Dict("label" => g.cands[jb].label, "k" => length(g.cands[jb].params), "bic" => g.cands[jb].bic, "isOptimal" => jb == 1)),
            "xmax" => g.kind == :s1b ? 20.0 : XMAX,
            "ymax" => g.kind == :s1b ? ceil(maximum(maximum(s.y) for s in series) * 1.2; sigdigits=2) : YMAX))
    end
    bundle = Dict{String,Any}("generated" => Dates.format(now(), "yyyy-mm-dd HH:MM"), "replicates" => B, "selection" => "fewest free parameters within 2 BIC of the lowest", "ruleBIC" => RULE_BIC, "maxCompartments" => MAX_COMPARTMENTS, "weighting" => selection_cfg("residual_weighting", "none"),
                              "zoo" => "autonomous", "groups" => out)
    mkpath(RESULTS_DIR)
    write(joinpath(RESULTS_DIR, "bundle.js"), stage_names("window.RESULTS = " * jsonv(bundle) * ";\n"))
    println("Wrote $(joinpath(RESULTS_DIR, "bundle.js"))")
end
