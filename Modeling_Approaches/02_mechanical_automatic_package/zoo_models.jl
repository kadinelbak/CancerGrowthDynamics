# Autonomous model zoo: no model uses time t directly.
#
# Every right-hand side depends only on the state and on two inputs: the cisplatin dose D (µM)
# and the low-resource level w_R (u_R = 0.955 at 0.5% FBS, 0 at 10% FBS). Each input of an
# environment is a constant or a piecewise-constant Schedule D(t) (data/schedules.toml). Delay and
# accumulation live in extra states driven by those inputs, so the same equations can be driven
# by any dose schedule D(t) — what optimal control needs. For a constant input switched on at
# t = 0 some of these reproduce the archived time-based forms exactly (an effect compartment
# gives D(1 − e^{−kt})), but they stay correct when the dose changes.
#
# One cell line in one environment is a "line": its states are
#   [X (observed count), damaged-cell compartments…, effect states…, auxiliary states…]
# growth   :logistic | :theta | :gompertz, using viable cells V = X − damaged, crowding by totals
# kernel   how the input w becomes the effect E
#   :none    E = w                                         immediate
#   :comp1   dA/dt = k (w − A),            E = A           effect compartment
#   :chain3  A1 → A2 → A3 at rate k,       E = A3          3 transit compartments (gamma delay)
#   :chain2  the same with 2 compartments  E = A2
#   :chain8  the same with 8 compartments  (supported, but not in the catalogue: > 3 compartments)
#   :second  d²A/dt² = ω²(w − A) − 2ζω dA/dt, E = max(A,0) second-order ("double derivative")
#   :auc     dA/dt = w,                    E = A           cumulative exposure, never cleared
# action   what the effect does to the cells
#   :loss      death (d0·w + β·E) per viable cell
#   :hill      death β E^h/(EC50^h + E^h), EC50 per line (the measured IC50, stages.toml [lines])
#   :damage    viable → damaged at β·E; damaged cells die at δ (still counted, do not divide)
#   :transit   as :damage, but damaged cells pass two compartments at rate τ before dying
#   :capacity  K_eff = K/(1 + c·E)
#   :fade      r_eff = r/(1 + c·E), plus immediate loss d0·w
#   :rho       r_eff = r(1 + (ρr − 1)w/u_R), K_eff = K(1 + (ρK − 1)w/u_R)
#   :sep_r     r_eff = r + (r_LR − r)w/u_R
#   :resource  R consumed at c·E·R·X/1000 (R(0) = 1); growth × R, death β(1 − R)
#   :toxin     W made at E·X/1000, cleared at κ; death γ·W

const SLOT_KEYS = (:r, :K, :θ, :k, :ω, :ζ, :β, :h, :δ, :d0, :c, :ρr, :ρK, :rLR, :γ, :κ, :τ, :ε, :d, :EC50)
# EC50 has no default: every line sets its own (a Hill line without one is an error, see build_system)
const SLOT_DEFAULT = (r=NaN, K=NaN, θ=1.0, k=1.0, ω=1.0, ζ=1.0, β=0.0, h=1.0, δ=1.0, d0=0.0, c=0.0, ρr=1.0, ρK=1.0,
                      rLR=NaN, γ=0.0, κ=0.0, τ=1.0, ε=1.0, d=0.0, EC50=NaN)

struct Slot
    i::Int        # index into the parameter vector, 0 = fixed
    v::Float64    # fixed value
end
@inline sv(s::Slot, p, ::Type{T}) where {T} = s.i > 0 ? T(p[s.i]) : T(s.v)
# EC50 of a module; lines saved before EC50 became a per-line slot (2026-10-02) used EC50 = 1 µM
@inline ec50_of(s, p, ::Type{T}) where {T} = hasproperty(s, :EC50) ? sv(s.EC50, p, T) : T(EC50)

# A piecewise-constant input: value v[i] from t[i] until t[i+1] (t[1] = 0). Constant inputs stay numbers.
struct Schedule
    t::Vector{Float64}
    v::Vector{Float64}
end
@inline input_at(x::Real, t) = x
@inline input_at(s::Schedule, t) = s.v[max(searchsortedlast(s.t, t), 1)]
breakpoints(x::Real) = Float64[]
breakpoints(s::Schedule) = s.t[2:end]
# times where an input of the system jumps (tstops for the ODE solver)
sys_breakpoints(sys) = sort(unique(reduce(vcat, [vcat(breakpoints(e.D), breakpoints(e.R)) for e in sys.envs]; init=Float64[])))

ncell(a::Symbol) = a == :damage ? 1 : a == :transit ? 2 : 0
chainlen(k::Symbol) = k == :chain2 ? 2 : k == :chain3 ? 3 : k == :chain8 ? 8 : 0
nker(k::Symbol) = chainlen(k) > 0 ? chainlen(k) : k == :second ? 2 : k in (:comp1, :auc) ? 1 : 0
naux(a::Symbol) = a in (:resource, :toxin) ? 1 : 0
# A line may carry extra mechanism modules (L.extra), e.g. a low-resource module next to its drug
# module in Stage 4. Each module keeps its own slots and its own block of hidden states, in order:
# [X, hidden states of module 1, hidden states of module 2, …]. Lines without `extra` are unchanged.
modules(L) = hasproperty(L, :extra) ? vcat(Any[L], L.extra) : Any[L]
mwidth(M) = ncell(M.action) + nker(M.kernel) + naux(M.action)
lwidth(L) = 1 + sum(mwidth(M) for M in modules(L))

# line spec: growth, kernel, action, input (:D or :R), slots (NamedTuple of Slot over SLOT_KEYS)
function line_spec(growth, kernel, action, input; fixed=Dict{Symbol,Float64}(), free=Dict{Symbol,Int}())
    vals = Dict{Symbol,Slot}()
    for k in SLOT_KEYS
        vals[k] = haskey(free, k) ? Slot(free[k], 0.0) : Slot(0, get(fixed, k, getfield(SLOT_DEFAULT, k)))
    end
    (growth=growth, kernel=kernel, action=action, input=input, s=NamedTuple{SLOT_KEYS}(Tuple(vals[k] for k in SLOT_KEYS)))
end

# env: inputs (D, R) and the lines present (indices into the system's line list)
# system: lines, envs, alpha[l][m] Slot (crowding of l by m), state layout.
# The observed counts come first, one per (env, line) in order, so every candidate of a group
# observes the same state indices; each line's hidden states follow in a block of their own:
#   xi[e][j] + 1 = observed count, hb[e][j] + 1 … = hidden states [damaged…, effect…, auxiliary…]
function build_system(lines, envs, alpha)
    for L in lines, M in modules(L)
        M.action == :hill && hasproperty(M.s, :EC50) && M.s.EC50.i == 0 && isnan(M.s.EC50.v) &&
            error("a Hill kill line has no EC50: set it from stages.toml [lines]")
    end
    xi = Vector{Vector{Int}}(); hb = Vector{Vector{Int}}()
    nobs = sum(length(e.lines) for e in envs); o = 0; h = nobs
    for e in envs
        v = Int[]; w = Int[]
        for l in e.lines
            push!(v, o); o += 1
            push!(w, h); h += lwidth(lines[l]) - 1
        end
        push!(xi, v); push!(hb, w)
    end
    (lines=lines, envs=envs, alpha=alpha, xi=xi, hb=hb, n=h)
end

function sys_u0(sys, X0)            # X0[e][j] = starting count of the j-th line in env e
    u = zeros(sys.n)
    for (e, env) in enumerate(sys.envs), (j, l) in enumerate(env.lines)
        L = sys.lines[l]
        u[sys.xi[e][j]+1] = X0[e][j]
        off = sys.hb[e][j]
        for M in modules(L)
            M.action == :resource && (u[off + ncell(M.action) + nker(M.kernel) + 1] = 1.0)
            off += mwidth(M)
        end
    end
    u
end
sys_obs(sys) = [sys.xi[e][j] + 1 for (e, env) in enumerate(sys.envs) for j in eachindex(env.lines)]
# states that are cell counts or concentrations and must stay ≥ 0 (the second-order velocity may be negative)
function sys_clamp(sys)
    m = trues(sys.n)
    for (e, env) in enumerate(sys.envs), (j, l) in enumerate(env.lines)
        off = sys.hb[e][j]
        for M in modules(sys.lines[l])
            M.kernel == :second && (m[off + ncell(M.action) + 2] = false)
            off += mwidth(M)
        end
    end
    m
end

@inline function growth_term(g, reff, V, x, θ)
    x = max(x, zero(x))
    g == :logistic && return reff * V * (1 - x)
    g == :theta && return reff * V * (1 - x^θ)
    return reff * V * log(1 / max(x, oftype(x, 1e-9)))       # gompertz
end

# One mechanism module: its delay kernel (hidden states from o+1) and its action. Returns the
# effective growth rate and capacity (given r, K), the extra per-viable-cell loss, and the outflow
# of dying damaged cells. Used by single-module lines (line_rhs!) and multi-module lines.
# ks: Norton–Simon factor for the damaged-cell actions, whose kill is the transfer β·E·V into the damaged
# compartment rather than a loss rate (1 = no coupling)
function module_effect!(du, u, o, M, p, w, X, V, r, K, ::Type{T}; ks=one(T)) where {T}
    s = M.s; z = zero(T)
    nc = ncell(M.action); nk = nker(M.kernel)
    w = w * sv(s.ε, p, T)
    ko = o + nc
    E = w
    if M.kernel == :comp1
        k = sv(s.k, p, T); du[ko+1] = k*(w - u[ko+1]); E = max(u[ko+1], z)
    elseif nk > 0 && chainlen(M.kernel) > 0
        k = sv(s.k, p, T)
        du[ko+1] = k*(w - u[ko+1])
        for j in 2:nk
            du[ko+j] = k*(u[ko+j-1] - u[ko+j])
        end
        E = max(u[ko+nk], z)
    elseif M.kernel == :second
        ω = sv(s.ω, p, T); ζ = sv(s.ζ, p, T)
        du[ko+1] = u[ko+2]; du[ko+2] = ω^2*(w - u[ko+1]) - 2ζ*ω*u[ko+2]
        E = max(u[ko+1], z)
    elseif M.kernel == :auc
        du[ko+1] = w; E = max(u[ko+1], z)
    end
    ao = ko + nk
    reff, Keff, loss, out = r, K, z, z
    a = M.action
    if a == :loss
        loss += sv(s.d0, p, T)*w + sv(s.β, p, T)*E
    elseif a == :hill
        h = sv(s.h, p, T)
        ec = ec50_of(s, p, T)
        loss += E > 0 ? sv(s.β, p, T) * E^h / (ec^h + E^h) : z
    elseif a == :damage
        δ = sv(s.δ, p, T); Dm = max(u[o+1], z)
        du[o+1] = ks*sv(s.β, p, T)*E*V - δ*Dm; out = δ*Dm
    elseif a == :transit
        τ = sv(s.τ, p, T); D1 = max(u[o+1], z); D2 = max(u[o+2], z)
        du[o+1] = ks*sv(s.β, p, T)*E*V - τ*D1; du[o+2] = τ*(D1 - D2); out = τ*D2
    elseif a == :capacity
        Keff = K / (1 + sv(s.c, p, T)*E)
    elseif a == :fade
        reff = r / (1 + sv(s.c, p, T)*E); loss += sv(s.d0, p, T)*w
    elseif a == :rho
        f = w / U_R; reff = r*(1 + (sv(s.ρr, p, T) - 1)*f); Keff = K*(1 + (sv(s.ρK, p, T) - 1)*f)
    elseif a == :sep_r
        reff = r + (sv(s.rLR, p, T) - r)*w/U_R
    elseif a == :resource
        R = max(u[ao+1], z)
        du[ao+1] = -sv(s.c, p, T)*E*R*X/1000; reff = r*R; loss += sv(s.β, p, T)*(1 - R)
    elseif a == :toxin
        W = max(u[ao+1], z)
        du[ao+1] = E*X/1000 - sv(s.κ, p, T)*W; loss += sv(s.γ, p, T)*W
    end
    reff, Keff, loss, out
end

# Norton–Simon (optional, L.ns): the drug kill is scaled by the line's growth factor, 1 − x for
# logistic and 1 − x^θ for theta-logistic growth (x = crowding / capacity, floored at 0), so crowded,
# slowly cycling cells die less. No extra parameter. Lines without `ns` are unchanged.
ns_on(L) = hasproperty(L, :ns) && L.ns
@inline function growth_factor(g, x, θ)
    x = max(x, zero(x))
    g == :logistic && return max(1 - x, zero(x))
    g == :theta && return max(1 - x^θ, zero(x))
    error("Norton–Simon coupling needs logistic or theta-logistic growth")
end

function line_rhs!(du, u, xi, o, L, p, w, crowd, ::Type{T}) where {T}
    s = L.s; z = zero(T)
    X = max(u[xi+1], z)
    cells = z
    for j in 1:ncell(L.action)
        cells += max(u[o+j], z)
    end
    V = max(X - cells, z)
    r = sv(s.r, p, T); K = sv(s.K, p, T); θ = sv(s.θ, p, T)
    # Norton–Simon: the kill is scaled by the growth factor; for the damaged-cell actions the transfer into the
    # damaged compartment is scaled (the module's own capacity is K, as no drug action changes K)
    ns = ns_on(L) && L.input == :D
    ks = ns && L.action in (:damage, :transit) ? growth_factor(L.growth, crowd / K, θ) : one(T)
    reff, Keff, loss, out = module_effect!(du, u, o, L, p, w, X, V, r, K, T; ks=ks)
    ns && (loss *= growth_factor(L.growth, crowd / Keff, θ))
    du[xi+1] = growth_term(L.growth, reff, V, crowd / Keff, θ) - (sv(s.d, p, T) + loss)*V - out
    nothing
end

# A line with several modules (L.extra). Growth-rate and capacity effects multiply, losses and
# outflows add. Drug-module losses are scaled by (1 + χ·u_R/0.955) when the line has an interaction
# slot L.chi: χ is the relative change of the cisplatin kill at 0.5% FBS (drug × low resource).
function line_rhs_multi!(du, u, xi, o, L, p, env, crowd, ::Type{T}) where {T}
    mods = modules(L); z = zero(T)
    X = max(u[xi+1], z)
    cells = z; off = o
    for M in mods
        for j in 1:ncell(M.action)
            cells += max(u[off+j], z)
        end
        off += mwidth(M)
    end
    V = max(X - cells, z)
    s = L.s; r = sv(s.r, p, T); K = sv(s.K, p, T); θ = sv(s.θ, p, T)
    reff, Keff, loss, out = r, K, sv(s.d, p, T), z
    χ = hasproperty(L, :chi) ? sv(L.chi, p, T) : z
    off = o
    if !ns_on(L)
        for M in mods
            w = M.input == :D ? T(env.D) : T(env.R)
            re, Ke, lo, ou = module_effect!(du, u, off, M, p, w, X, V, r, K, T)
            reff *= re / r; Keff *= Ke / K
            loss += M.input == :D ? lo * (1 + χ * T(env.R) / T(U_R)) : lo
            out += ou; off += mwidth(M)
        end
    else                                     # Norton–Simon: drug losses scaled by the growth factor
        any(M -> M.input == :D && M.action in (:damage, :transit), mods) &&
            error("Norton–Simon with a damaged-cell drug action is only implemented for single-module lines")
        dloss = z
        for M in mods
            w = M.input == :D ? T(env.D) : T(env.R)
            re, Ke, lo, ou = module_effect!(du, u, off, M, p, w, X, V, r, K, T)
            reff *= re / r; Keff *= Ke / K
            M.input == :D ? (dloss += lo * (1 + χ * T(env.R) / T(U_R))) : (loss += lo)
            out += ou; off += mwidth(M)
        end
        loss += dloss * growth_factor(L.growth, crowd / Keff, θ)
    end
    du[xi+1] = growth_term(L.growth, reff, V, crowd / Keff, θ) - loss*V - out
    nothing
end

function make_system_model(sys)
    function f!(du, u, p, t)
        T = promote_type(eltype(u), eltype(p))
        for (e, envs) in enumerate(sys.envs)
            env = (D=input_at(envs.D, t), R=input_at(envs.R, t), lines=envs.lines)      # inputs at time t
            for (j, l) in enumerate(env.lines)
                xi = sys.xi[e][j]; L = sys.lines[l]
                crowd = max(u[xi+1], zero(T))
                for (m, lm) in enumerate(env.lines)
                    m == j && continue
                    crowd += sv(sys.alpha[l][lm], p, T) * max(u[sys.xi[e][m]+1], zero(T))
                end
                if hasproperty(L, :extra)
                    line_rhs_multi!(du, u, xi, sys.hb[e][j], L, p, env, crowd, T)
                else
                    w = L.input == :D ? T(env.D) : T(env.R)
                    line_rhs!(du, u, xi, sys.hb[e][j], L, p, w, crowd, T)
                end
            end
        end
        nothing
    end
end

# ── symbolic equations ─────────────────────────────────────────────────────────

const GREEK = Dict("θ" => "\\theta", "ρ" => "\\rho", "β" => "\\beta", "λ" => "\\lambda", "α" => "\\alpha", "δ" => "\\delta",
                   "γ" => "\\gamma", "κ" => "\\kappa", "τ" => "\\tau", "ω" => "\\omega", "ζ" => "\\zeta", "ε" => "\\varepsilon",
                   "χ" => "\\chi")
function ptex(n)
    base, rest = n, ""
    if occursin("^", n)
        base, sup = split(n, "^"; limit=2); rest = "^{$sup}"
    end
    parts = split(base, "_"; limit=2)
    b = get(GREEK, string(parts[1]), string(parts[1]))
    length(parts) == 2 ? "$(b)_{$(parts[2])}$rest" : "$b$rest"
end

# display name of a slot for line symbol X and input (:D drug, :R low resource)
function slot_name(k::Symbol, X, input; combo=false)
    if input == :D
        return Dict(:β => "β_$X", :k => "k_$X", :h => "h_$X", :δ => "δ_$X", :τ => "τ_$X", :ω => "ω_$X", :ζ => "ζ_$X",
                    :ε => "ε_$X", :d => "d_$X", :r => "r_$X", :K => "K_$X", :θ => "θ_$X", :d0 => "d_$(X)^0")[k]
    end
    Dict(:β => combo ? "d_$(X)^acc" : "d_$(X)^LR", :d0 => "d_$(X)^LR", :k => "k_$(X)^R", :h => "h_$(X)^R", :δ => "δ_$(X)^R",
         :τ => "τ_$(X)^R", :ω => "ω_$(X)^R", :ζ => "ζ_$(X)^R", :c => "c_$X", :ρr => "ρ_r", :ρK => "ρ_K", :rLR => "r_$(X)^LR",
         :γ => "γ_$X", :κ => "κ_$X", :ε => "ε_$X", :d => "d_$X", :r => "r_$X", :K => "K_$X", :θ => "θ_$X")[k]
end

# TeX lines for one line. X, Y: symbols of this line and the other line (Y = "" in monoculture);
# crowd: TeX of the crowding numerator; Ksym: TeX of K
isoff(s::Slot) = s.i == 0 && s.v == 0
# TeX pieces of one module: its delay-state equations, how it changes r and K, its loss terms
# (each a string starting with "-"), and extra state equations. V: symbol of viable cells.
function module_tex(M, X, V)
    # a mechanism whose input or strength is fixed at 0 is written as plain growth
    if isoff(M.s.ε) || (M.action in (:loss, :hill, :damage, :transit) && isoff(M.s.β) && isoff(M.s.d0))
        M = merge(M, (kernel=:none, action=:none))
    end
    input = M.input; combo = M.s.d0.i > 0 || M.s.d0.v != 0
    nm(k) = ptex(slot_name(k, X, input; combo=combo))
    w = input == :D ? "D" : "u_R"
    M.s.ε.i > 0 && (w = "$(nm(:ε))\\,$w")
    Asym = input == :D ? "A_{$X}" : "S_{$X}"
    lines = String[]; losses = String[]
    E = w
    if M.kernel == :comp1
        push!(lines, "\\frac{d$Asym}{dt}=$(nm(:k))\\left($w-$Asym\\right)"); E = Asym
    elseif chainlen(M.kernel) > 0
        n = chainlen(M.kernel)
        a(i) = input == :D ? "A_{$X,$i}" : "S_{$X,$i}"
        push!(lines, "\\frac{d$(a(1))}{dt}=$(nm(:k))\\left($w-$(a(1))\\right),\\ \\frac{d$(a("j"))}{dt}=$(nm(:k))\\left($(a("j-1"))-$(a("j"))\\right),\\ j=2,\\dots,$n")
        E = a(n)
    elseif M.kernel == :second
        push!(lines, "\\frac{d^2$Asym}{dt^2}={$(nm(:ω))}^2\\left($w-$Asym\\right)-2$(nm(:ζ))$(nm(:ω))\\frac{d$Asym}{dt}"); E = Asym
    elseif M.kernel == :auc
        push!(lines, "\\frac{d$Asym}{dt}=$w"); E = Asym
    end
    a = M.action
    rfun(r) = a == :fade ? "\\frac{$r}{1+$(nm(:c))$E}" : a == :rho ? "\\rho_r $r" : a == :sep_r ? nm(:rLR) : a == :resource ? "$r R_{$X}" : r
    Kfun(K) = a == :capacity ? "\\frac{$K}{1+$(nm(:c))$E}" : a == :rho ? "\\rho_K $K" : K
    if a == :loss
        push!(losses, M.s.d0.i > 0 || M.s.d0.v != 0 ? "-\\left($(nm(:d0))$w+$(nm(:β))$E\\right)$V" : "-$(nm(:β))\\,$E\\,$V")
    elseif a == :hill
        push!(losses, "-$(nm(:β))\\frac{$E^{$(nm(:h))}}{EC_{50,$X}^{$(nm(:h))}+$E^{$(nm(:h))}}\\,$V")
    elseif a == :damage
        push!(losses, "-$(nm(:δ))\\,$(X)_D")
        push!(lines, "\\frac{d$(X)_D}{dt}=$(nm(:β))\\,$E\\,$V-$(nm(:δ))\\,$(X)_D,\\quad $V=$X-$(X)_D")
    elseif a == :transit
        push!(losses, "-$(nm(:τ))\\,$(X)_{D2}")
        push!(lines, "\\frac{d$(X)_{D1}}{dt}=$(nm(:β))\\,$E\\,$V-$(nm(:τ))$(X)_{D1},\\ \\frac{d$(X)_{D2}}{dt}=$(nm(:τ))\\left($(X)_{D1}-$(X)_{D2}\\right),\\quad $V=$X-$(X)_{D1}-$(X)_{D2}")
    elseif a == :fade
        (M.s.d0.i > 0 || M.s.d0.v != 0) && push!(losses, "-$(nm(:d0))\\,$w\\,$X")
    elseif a == :resource
        push!(losses, "-$(nm(:β))\\left(1-R_{$X}\\right)$X")
        push!(lines, "\\frac{dR_{$X}}{dt}=-$(nm(:c))\\,$E\\,R_{$X}\\frac{$X}{1000},\\ R_{$X}(0)=1")
    elseif a == :toxin
        push!(losses, "-$(nm(:γ))\\,W_{$X}\\,$X")
        push!(lines, "\\frac{dW_{$X}}{dt}=$E\\frac{$X}{1000}" * (M.s.κ.i > 0 || M.s.κ.v != 0 ? "-$(nm(:κ))W_{$X}" : ""))
    end
    (lines=lines, losses=losses, rfun=rfun, Kfun=Kfun, input=input)
end

# TeX lines for one line (any number of modules). X, Y: symbols of this line and the other line;
# crowd: TeX of the crowding numerator; Ksym: TeX of K
function line_tex(L, X, crowd, Ksym)
    mods = modules(L)
    V = any(M -> M.action in (:damage, :transit), mods) ? "V_{$X}" : X
    parts = [module_tex(M, X, V) for M in mods]
    r = ptex(slot_name(:r, X, L.input)); K = Ksym; th = ptex(slot_name(:θ, X, L.input))
    for q in parts
        r = q.rfun(r); K = q.Kfun(K)
    end
    frac = "\\frac{$crowd}{$K}"
    G = L.growth == :logistic ? "$r\\,$V\\left(1-$frac\\right)" :
        L.growth == :theta ? "$r\\,$V\\left[1-\\left($frac\\right)^{$th}\\right]" : "$r\\,$V\\,\\ln\\left(\\frac{$K}{$crowd}\\right)"
    rhs = G
    (L.s.d.i > 0 || L.s.d.v != 0) && (rhs *= "-$(ptex(slot_name(:d, X, L.input)))\\,$V")
    chi = hasproperty(L, :chi) && !isoff(L.chi)
    nsf = !ns_on(L) ? "" : L.growth == :theta ? "\\left[1-\\left($frac\\right)^{$th}\\right]" : "\\left(1-$frac\\right)"
    for q in parts, lo in q.losses
        lo = (chi && q.input == :D) ? "-\\left(1+\\chi_{$X}\\,\\hat u_R\\right)\\left[$(lo[2:end])\\right]" : lo
        rhs *= (nsf != "" && q.input == :D) ? "-$nsf\\left[$(lo[2:end])\\right]" : lo
    end
    lines = reduce(vcat, [q.lines for q in parts]; init=String[])
    chi && push!(lines, "\\hat u_R=u_R/0.955\\ (1\\text{ at }0.5\\%\\text{ FBS})")
    vcat(["\\frac{d$X}{dt}=$rhs"], lines)
end

# ── candidate catalogue ────────────────────────────────────────────────────────
# (name, label, family, kernel, action, free slots => (p0 list, bounds))

# The catalogue is a grid: every delay kernel × every effect, plus a few special mechanisms.
const KERNELS = [(:none, "imm", "Immediate"), (:comp1, "comp", "Effect compartment"), (:chain2, "chain2", "2 transit compartments"),
                 (:chain3, "chain3", "3 transit compartments"),
                 (:second, "second", "Second order (double derivative)"), (:auc, "auc", "Cumulative exposure (no clearance)")]
# At most 3 compartments per model: delay states plus damaged-cell states (the second-order effect counts as 2).
const MAX_COMPARTMENTS = 3
# selection: fewest free parameters among the candidates within RULE_BIC of the lowest BIC
const RULE_BIC = 2.0   # user, 2026-10-27 (was 3.0)
ncompartments(kernel, action) = nker(kernel) + ncell(action) + naux(action)
kernel_slots(k) = chainlen(k) > 0 || k == :comp1 ? [:k] : k == :second ? [:ω, :ζ] : Symbol[]
action_slots(a) = Dict(:loss => [:β], :hill => [:β, :h], :damage => [:β, :δ], :transit => [:β, :τ], :capacity => [:c],
                       :fade => [:c, :d0], :rho => [:ρr, :ρK], :sep_r => [:rLR], :resource => [:c, :β], :toxin => [:γ])[a]
const KERNEL_WORD = Dict(:none => "", :comp1 => "effect compartment", :chain2 => "2 compartments", :chain3 => "3 compartments",
                         :chain8 => "8 compartments", :second => "second-order effect", :auc => "cumulative exposure")
function zoo_grid(input)
    word = input == :D ? Dict(:loss => "linear kill", :hill => "Hill kill", :damage => "damaged cells die", :capacity => "lowers K") :
                         Dict(:loss => "loss", :hill => "Hill loss", :damage => "damaged cells die", :capacity => "lowers K")
    acts = input == :D ? [:loss, :hill, :damage] : [:loss, :damage, :capacity]
    pre = input == :D ? "" : "lr_"
    out = Tuple{String,String,String,Symbol,Symbol,Vector{Symbol}}[]
    for (k, short, fam) in KERNELS, a in acts
        ncompartments(k, a) <= MAX_COMPARTMENTS || continue
        label = k == :none ? "Immediate " * word[a] : "$(uppercasefirst(KERNEL_WORD[k])) → $(word[a])"
        push!(out, ("$(pre)$(short)_$(a)", label, fam, k, a, vcat(action_slots(a), kernel_slots(k))))
    end
    # immediate + delayed loss together
    for k in (:comp1, :chain2, :chain3)
        push!(out, ("$(pre)imm_plus_$(k)_loss", "Immediate + $(KERNEL_WORD[k]) $(word[:loss])", "Immediate + delayed", k, :loss, [:d0, :β, :k]))
    end
    out
end

const DRUG_ZOO = zoo_grid(:D)

const LR_ZOO = vcat(zoo_grid(:R), [
    ("lr_rho",          "Growth-rate and K multipliers",                   "Rescaled growth",       :none,  :rho,      [:ρr, :ρK]),
    ("lr_sep_r",        "Separate low-resource growth rate",               "Rescaled growth",       :none,  :sep_r,    [:rLR]),
    ("lr_comp_fade",    "Effect compartment slows growth, immediate loss", "Effect compartment",    :comp1, :fade,     [:c, :k, :d0]),
    ("lr_chain2_fade",  "2 compartments slow growth, immediate loss",      "2 transit compartments",:chain2,:fade,     [:c, :k, :d0]),
    ("lr_transit",      "Damaged cells die through two stages",            "Damaged cells",         :none,  :transit,  [:β, :τ]),
    ("lr_comp_transit", "Effect compartment → two-stage damage → death",   "Damaged cells",         :comp1, :transit,  [:β, :k, :τ]),
    ("lr_resource",     "Cells use up the resource, then starve",          "Resource / by-product", :none,  :resource, [:c, :β]),
    ("lr_toxin",        "Toxic by-product accumulates",                    "Resource / by-product", :none,  :toxin,    [:γ]),
    ("lr_toxin_clear",  "Toxic by-product accumulates and clears",         "Resource / by-product", :none,  :toxin,    [:γ, :κ]),
])
@assert all(ncompartments(c[4], c[5]) <= MAX_COMPARTMENTS for c in vcat(DRUG_ZOO, LR_ZOO))

# ζ ≥ 1 (user, 2026-11-02, P5 of the decision review): the second-order effect is critically damped or overdamped,
# so it rises and settles without oscillating (ζ_N sat on its old lower bound of 0.05, an effect that rang like an
# undamped spring to bend through one dip in the A2780 low-serum counts)
const SLOT_BOUNDS = Dict(:β => (0.0, 50.0), :k => (1e-3, 20.0), :h => (0.1, 10.0), :δ => (1e-3, 10.0), :τ => (1e-3, 10.0),
                         :ω => (1e-3, 20.0), :ζ => (1.0, 5.0), :d0 => (0.0, 50.0), :c => (0.0, 1000.0), :ρr => (0.0, 20.0),
                         :ρK => (1e-3, 3.0), :rLR => (1e-4, 5.0), :γ => (0.0, 10.0), :κ => (0.0, 10.0), :ε => (0.0, 5.0), :d => (0.0, 3.0))
const SLOT_STARTS = Dict(:β => [0.3, 1.0, 0.05], :k => [0.5, 2.0, 0.1], :h => [2.0, 1.0, 4.0], :δ => [0.5, 2.0, 0.1], :τ => [0.5, 2.0, 0.1],
                         :ω => [0.5, 2.0, 0.2], :ζ => [1.0, 2.0, 4.0], :d0 => [0.3, 0.05, 1.0], :c => [1.0, 5.0, 0.2], :ρr => [1.0, 2.0, 0.5],
                         :ρK => [0.8, 0.3, 1.2], :rLR => [0.5, 0.1, 1.0], :γ => [0.1, 0.5, 0.02], :κ => [0.1, 0.5, 0.01], :ε => [1.0, 0.5, 2.0], :d => [0.1, 0.02, 0.5])
