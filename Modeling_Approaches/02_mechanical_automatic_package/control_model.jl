# Control plant for optimal control of cisplatin dosing, built from the selected compartment models.
#
# The biology is NOT refitted here. The cell equations are the zoo's own right-hand side
# (zoo_models.jl: line_rhs!; simulated through a flattened, type-stable copy that the tests check against
# it to 1e-12) with the parameters of the carried-forward models in outputs/zoo_fits/zoo_state.jls.
# Only the dose input and a toxicity proxy are added:
#
#   state x = [N, C, A_N,1..3, A_C,1..3, (low resource: S_C,1..3), B, T]
#     N, C      sensitive (A2780) and resistant (A2780cis) cells (Mean Cells units)
#     A_X,j     drug-effect compartments of line X (3-compartment chain, rate k_X, driven by u)
#     S_C,j     low-resource stress compartments of A2780cis (low-resource scenario only)
#     B         systemic drug burden, dB/dt = s_B·u − d_B·B               (de Pillis & Radunskaya 2001)
#     T         toxicity burden 1 − H, dT/dt = a_H(1 − e^(−B))(1 − T) − r_H·T(1 − T)
#   control u = cisplatin in the medium (µM), piecewise constant over each decision interval
#
# Scenarios for the drug response in coculture (user decision, 2026-09-29):
#   :A  the Stage 3 optimal model as selected (coculture fit; the resistant kill is inherited from
#       Stage 2, refitted, or fixed at 0 when the selected candidate has no resistant cisplatin term:
#       stage3_betaC() says which)
#   :B  Stage 2 monoculture drug response for both lines, crowding fitted to untreated cocultures
#       only (the variant that predicted the held-out treated cocultures best)

include(joinpath(@__DIR__, "zoo_common.jl"))
include(joinpath(@__DIR__, "zoo_models.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization, TOML, ForwardDiff, Statistics
import LinearAlgebra
const LA = LinearAlgebra

const CFG = TOML.parsefile(joinpath(@__DIR__, "control_config.toml"))
# ZOO_STATE / AN_STATE override the fitted states (e.g. a backed-up earlier fit, to compare control runs)
const ZOO_STATE = deserialize(get(ENV, "ZOO_STATE", joinpath(@__DIR__, "outputs", "zoo_fits", "zoo_state.jls")))
const AN_STATE = deserialize(get(ENV, "AN_STATE", joinpath(@__DIR__, "outputs", "zoo_analyses", "analyses_state.jls")))
zgroup(id) = ZOO_STATE.groups[findfirst(g -> g.id == id, ZOO_STATE.groups)]

# ── reuse of finished control runs (2026-11-01) ──────────────────────────────────────────────────────────
# A control step (run_control.jl, run_control_fasting.jl, run_control_sweep.jl) is skipped when everything it
# reads is unchanged since its last finished run: the control code and configs, the data, the package
# versions, and the parts of the fitted states it uses (each group's carried model and its bootstrap draws,
# the Ce0 crowding fit). A refit that leaves the carried models and draws the same therefore costs no control
# time; anything else reruns. CONTROL_FORCE=1 (run_pipeline.jl --recompute-control) always recomputes.
# The fingerprint is a SHA-256 written next to the outputs after a successful run.
const CONTROL_CODE = ["control_model.jl", "control_mpc.jl", "control_fasting.jl", "zoo_common.jl", "zoo_models.jl",
                      "zoo_data.jl", "bootstrap_tools.jl", "stage_names.jl", "control_config.toml", "stages.toml",
                      "Manifest.toml", joinpath("data", "schedules.toml"), joinpath("data", "measurements.csv")]
function control_fingerprint(script::AbstractString, extra::Pair...)
    io = IOBuffer()
    for f in vcat(CONTROL_CODE, basename(script))
        fp = joinpath(@__DIR__, f)
        write(io, f, 0x00, isfile(fp) ? read(fp) : UInt8[], 0x00)
    end
    for g in ZOO_STATE.groups
        c = g.cands[1]
        write(io, string(g.id), 0x00, string(c.name), 0x00, join(string.(c.pnames), ","), 0x00, Float64.(c.params))
        if hasproperty(c, :fixedcc)     # sorted: a Dict's order is not stable between sessions
            fc = c.fixedcc isa AbstractDict ? sort!([string(k) => v for (k, v) in c.fixedcc]; by=first) : c.fixedcc
            write(io, repr(fc), 0x00)
        end
        if haskey(ZOO_STATE.results, g.id)
            P = ZOO_STATE.results[g.id].res.pars[1]
            write(io, Int64(size(P, 1)), Int64(size(P, 2)), Float64.(vec(P)))
        end
    end
    write(io, repr(ZOO_STATE.B), 0x00, Float64.(AN_STATE[:validation].ce1[1].params))
    for (k, v) in extra
        write(io, string(k), "=", string(v), 0x00)
    end
    bytes2hex(sha256(take!(io)))
end
# true when `fpfile` holds this fingerprint and every output in `outs` exists (and CONTROL_FORCE is not 1)
control_unchanged(fpfile, fp, outs...) = get(ENV, "CONTROL_FORCE", "0") != "1" && isfile(fpfile) &&
    strip(read(fpfile, String)) == fp && all(isfile, outs)

# resistant kill β_C of the carried Stage 3 model, and where it comes from: refitted, dropped (the
# "no resistant cisplatin term" candidates fix it at 0) or inherited from Stage 2b
function stage3_betaC()
    c3 = zgroup("s3").cands[1]; s = c3.sys.lines[2].s.β
    dropped = s.i == 0 && hasproperty(c3, :fixedcc) && haskey(c3.fixedcc, :βC)
    s.i > 0 ? (value=c3.params[s.i], free=true, dropped=false,
               source="Stage 3 (refitted" * (c3.bounds[s.i][1] > 0 ? ", floor $(round(c3.bounds[s.i][1]; sigdigits=3)))" : ", no floor)")) :
    dropped ? (value=s.v, free=false, dropped=true, source="Stage 3 (no resistant cisplatin term in coculture: selected, β_C = 0)") :
              (value=s.v, free=false, dropped=false, source="Stage 2 (inherited)")
end

# ── flattened cell equations ─────────────────────────────────────────────────
# Flattened, type-stable copy of the plant's cell equations (the same modules and values as the zoo's
# line_rhs! / line_rhs_multi!, which are generic and dispatch dynamically, about 50× slower). The tests
# check it against the zoo right-hand side to 1e-12 (normal medium: test_control.jl; with the low-resource
# modules: test_control_fasting.jl). Only the module types of the carried models are supported;
# line_specs stops with an error on anything else rather than silently differing.
struct ModSpec
    input_D::Bool; kernel_chain::Int; k::Float64; o::Int     # chain length (comp1 = 1, chain2 = 2, chain3 = 3; 0 = none/auc/second), rate, first hidden index − 1
    action::Symbol; β::Float64; h::Float64; c::Float64; d0::Float64; ec50::Float64
    ε::Float64; εj::Int                                       # fixed strength, or index into the free strengths (0 = fixed)
    kernel::Symbol; ω::Float64; ζ::Float64                    # delay kernel (zoo_models.jl KERNELS) and second-order ω, ζ
end
struct LineSpec
    ix::Int; other::Int; α::Float64; logistic::Bool; r::Float64; K::Float64; θ::Float64; d::Float64
    mods::Vector{ModSpec}
    ns::Bool                    # Norton–Simon: drug losses × growth factor (zoo_models.jl growth_factor)
    chi::Float64                # drug × low-resource interaction χ: drug losses × (1 + χ·u_R/U_R) (zoo line_rhs_multi!)
end
const FAST_KERNELS = (:none, :comp1, :chain2, :chain3, :second, :auc)
# iε: positions in p of free low-resource strengths (control_fasting.jl); none in normal medium
function line_specs(S, p, iε)
    out = LineSpec[]
    for (j, l) in enumerate(S.envs[1].lines)
        L = S.lines[l]; xi = S.xi[1][j]; off = S.hb[1][j]
        L.growth in (:logistic, :theta) || error("growth law $(L.growth) not supported in the fast plant")
        mods = ModSpec[]
        for M in modules(L)
            ncell(M.action) == 0 && naux(M.action) == 0 || error("action $(M.action) not supported in the fast plant")
            M.kernel in FAST_KERNELS || error("kernel $(M.kernel) not supported in the fast plant")
            M.action in (:hill, :capacity, :loss) || error("action $(M.action) not supported in the fast plant")
            nk = nker(M.kernel); gs(k) = hasproperty(M.s, k) ? sv(getfield(M.s, k), p, Float64) : 0.0
            εs = M.s.ε; εj = εs.i > 0 ? findfirst(==(εs.i), iε) : 0
            εs.i > 0 && εj === nothing && error("unexpected free slot")
            ch = M.kernel in (:comp1, :chain2, :chain3) ? nk : 0
            push!(mods, ModSpec(M.input == :D, ch, ch > 0 ? gs(:k) : 0.0, off, M.action, gs(:β), gs(:h), gs(:c), gs(:d0), ec50_of(M.s, p, Float64),
                                εs.i > 0 ? NaN : εs.v, εj, M.kernel, gs(:ω), gs(:ζ)))
            off += mwidth(M)
        end
        m = 3 - j
        push!(out, LineSpec(xi + 1, S.xi[1][m] + 1, sv(S.alpha[l][S.envs[1].lines[m]], p, Float64), L.growth == :logistic,
                            sv(L.s.r, p, Float64), sv(L.s.K, p, Float64), sv(L.s.θ, p, Float64), sv(L.s.d, p, Float64), mods, ns_on(L),
                            hasproperty(L, :chi) ? sv(L.chi, p, Float64) : 0.0))
    end
    out
end
# Cell (and drug/stress chain) derivatives for dose u and low-resource level uR. f = ε factor of the dose,
# εv = free low-resource strengths, koff/v = refeeding fade rate (NaN = fitted rate) and fasting input.
function spec_cells_rhs!(du, x, u, uR, f, spec, εv, koff, v)
    T = promote_type(eltype(x), typeof(u), typeof(uR)); z = zero(T)
    @inbounds for L in spec
        X = max(x[L.ix], z); crowd = X + L.α * max(x[L.other], z)
        Keff = T(L.K); loss = T(L.d); dloss = z                   # dloss: drug-module losses (Norton–Simon scales these)
        for M in L.mods
            ε = M.εj == 0 ? T(M.ε) : 1 + (εv[M.εj] - 1) * f
            wi = (M.input_D ? T(u) : T(uR)) * ε
            E = wi; o = M.o
            if M.kernel_chain > 0                    # effect compartment / transit chain of n compartments at rate k
                # low-resource stress may fade faster after refeeding (sensitivity only; default = fitted rate)
                k = (M.input_D || isnan(koff)) ? T(M.k) : M.k * v + koff * (1 - v)
                du[o+1] = k * (wi - x[o+1])
                for i in 2:M.kernel_chain
                    du[o+i] = k * (x[o+i-1] - x[o+i])
                end
                E = max(x[o+M.kernel_chain], z)
            elseif M.kernel == :auc                  # cumulative exposure: never fades in the fitted model; the
                fade = (M.input_D || isnan(koff)) ? z : koff * (1 - v)       # refeeding sensitivity adds a fade when fed
                du[o+1] = wi - fade * x[o+1]
                E = max(x[o+1], z)
            elseif M.kernel == :second               # second-order (damped oscillator) delay
                du[o+1] = x[o+2]; du[o+2] = M.ω^2 * (wi - x[o+1]) - 2M.ζ * M.ω * x[o+2]
                E = max(x[o+1], z)
            end
            if M.action == :hill
                lo = E > 0 ? M.β * E^M.h / (T(M.ec50)^M.h + E^M.h) : z
                M.input_D ? (dloss += lo) : (loss += lo)
            elseif M.action == :capacity
                Keff /= 1 + M.c * E
            else                                     # :loss
                lo = M.d0 * wi + M.β * E
                M.input_D ? (dloss += lo) : (loss += lo)
            end
        end
        q = max(crowd / Keff, z)
        dloss *= 1 + L.chi * T(uR) / T(U_R)        # drug × low-resource interaction (0 unless Stage 4 carries χ)
        loss += L.ns ? dloss * (L.logistic ? max(1 - q, z) : max(1 - q^L.θ, z)) : dloss
        du[L.ix] = L.r * X * (L.logistic ? 1 - q : 1 - q^L.θ) - loss * X
    end
    nothing
end
const NO_EPS = Float64[]

# ── plant ─────────────────────────────────────────────────────────────────────

struct Plant
    name::String
    sys::Any                  # zoo system (one environment, lines [N, C])
    p::Vector{Float64}        # free-slot parameter vector of the zoo lines
    nx::Int
    iN::Int; iC::Int; iB::Int; iT::Int
    K::Float64                # shared carrying capacity (for normalisation and the cancer bound)
    tox::NamedTuple           # r_H, a_H, d_B, s_B
    uR::Float64               # low-resource level of the environment (0 in normal medium)
    fp::Any                   # low-resource environment: the two-input plant of control_fasting.jl held at v = 1, else nothing
    spec::Any                 # normal medium: flattened cell equations (spec_cells_rhs!), or nothing = the zoo right-hand side
end
line_specs(base::Plant, iε) = line_specs(base.sys, base.p, iε)

# replace fixed slot values of a line (e.g. a bootstrap draw of carried parameters)
function set_slots(L, d::Dict{Symbol,Float64})
    s = L.s
    for (k, v) in d
        s = merge(s, NamedTuple{(k,)}((Slot(0, v),)))
    end
    merge(L, (s=s,))
end

default_tox() = (r_H=Float64(CFG["toxicity"]["r_H"]), a_H=Float64(CFG["toxicity"]["a_H"]),
                 d_B=Float64(CFG["toxicity"]["d_B"]), s_B=Float64(CFG["toxicity"]["s_B"]))

# fast = true: simulate the cells with the flattened equations (normal medium only; a model with modules the
# flattened copy does not support falls back to the zoo right-hand side, which gives the same values, slower)
function make_plant(name, lines, alpha, p; uR=0.0, tox=default_tox(), fast=true)
    sys = build_system(lines, [(D=0.0, R=uR, lines=[1, 2])], alpha)
    n = sys.n
    Kslot = lines[1].s.K
    K = Kslot.i > 0 ? p[Kslot.i] : Kslot.v
    spec = (fast && uR == 0) ? (try line_specs(sys, p, Int[]) catch; nothing end) : nothing
    Plant(name, sys, p, n + 2, sys.xi[1][1] + 1, sys.xi[1][2] + 1, n + 1, n + 2, K, tox, uR, nothing, spec)
end

# The Stage 3 coculture lines and parameters for a scenario (point estimates or a bootstrap draw b)
function scenario_plant(scen::Symbol; b::Union{Nothing,Int}=nothing, tox=default_tox())
    g3 = zgroup("s3"); c3 = g3.cands[1]
    LN, LC = c3.sys.lines
    # a carried parameter of an earlier stage by name ("θ" → θ_N or θ_C), from its fit or bootstrap replicate b;
    # nothing when the carried model has no such parameter (e.g. h for a linear kill)
    function draw(gid, sym)
        c = zgroup(gid).cands[1]; i = findfirst(n -> n in ("$(sym)_N", "$(sym)_C"), c.pnames)
        i === nothing && return nothing
        b === nothing ? c.params[i] : ZOO_STATE.results[gid].res.pars[1][b, i]
    end
    setdraw!(d, sym, gid) = (v = draw(gid, sym); v === nothing || (d[sym] = v); d)
    fixed_N = Dict{Symbol,Float64}(); fixed_C = Dict{Symbol,Float64}()
    if b !== nothing       # carried parameters from the same bootstrap replicate of their own stage
        setdraw!(fixed_N, :r, "s1-a2780"); setdraw!(fixed_N, :θ, "s1-a2780"); setdraw!(fixed_C, :r, "s1-a2780cis"); setdraw!(fixed_C, :θ, "s1-a2780cis")
        for s in (:h, :k)
            setdraw!(fixed_N, s, "s2-a2780"); setdraw!(fixed_C, s, "s2-a2780cis")
        end
    end
    if scen == :A
        p = b === nothing ? copy(c3.params) : ZOO_STATE.results["s3"].res.pars[1][b, :]
        # an inherited resistant kill follows its own Stage 2 bootstrap draw, as h and k do
        b !== nothing && LC.s.β.i == 0 && !stage3_betaC().dropped && setdraw!(fixed_C, :β, "s2-a2780cis")
        name = "A: Stage 3 coculture model as selected"
    else
        # crowding (K, α_NC, α_CN and d_N where the model has it) fitted to the untreated cocultures; these
        # lead the Stage 3 parameter vector (run_heldout.jl checks the same), the kills come from Stage 2.
        # The untreated-coculture crowding fit has no bootstrap of its own: a draw b BORROWS the Stage 3a
        # bootstrap ratio of replicate b to its estimate (pars[b] / estimate), applied to the Ce0 values
        # (an approximation, disclosed in the report)
        ce0 = AN_STATE[:validation].ce1[1].params; n = length(ce0)
        p = b === nothing ? copy(ce0) : ZOO_STATE.results["s3"].res.pars[1][b, 1:n] .* (ce0 ./ c3.params[1:n])
        setdraw!(fixed_N, :β, "s2-a2780"); setdraw!(fixed_C, :β, "s2-a2780cis")
        name = "B: monoculture drug response for both lines, crowding from untreated cocultures"
    end
    any(isnan, p) && return nothing
    make_plant(name, [set_slots(LN, fixed_N), set_slots(LC, fixed_C)], c3.sys.alpha, p; tox=tox)
end

# Low-resource environment: lr_plant(scen) is defined in control_fasting.jl (the two-input plant held at
# v = 1), so the 14-day and the long (dose + fasting) controllers use the same dose-dependent drug ×
# low-resource model.

# continuous-time right-hand side for a given dose u (µM); ref = true forces the zoo's own right-hand side
# (the reference the flattened equations are tested against)
function plant_rhs!(du, x, u, P::Plant; ref=false)
    T = promote_type(eltype(x), typeof(u))
    P.fp === nothing || return fast_rhs!(du, x, (T(u), one(T)), P.fp)      # low resource: fasting plant at v = 1
    if P.spec !== nothing && !ref                                          # normal medium, flattened equations
        spec_cells_rhs!(du, x, T(u), zero(T), zero(T), P.spec, NO_EPS, NaN, 0)
    else
        zoo_cells_rhs!(du, x, u, P, T)
    end
    B = x[P.iB]; Tx = x[P.iT]; tp = P.tox
    du[P.iB] = tp.s_B * u - tp.d_B * B
    du[P.iT] = tp.a_H * (1 - exp(-max(B, zero(T)))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
    nothing
end
function zoo_cells_rhs!(du, x, u, P::Plant, T)
    S = P.sys; env = S.envs[1]
    for (j, l) in enumerate(env.lines)
        xi = S.xi[1][j]; L = S.lines[l]
        crowd = max(x[xi+1], zero(T))
        for (m, lm) in enumerate(env.lines)
            m == j && continue
            crowd += sv(S.alpha[l][lm], P.p, T) * max(x[S.xi[1][m]+1], zero(T))
        end
        if hasproperty(L, :extra)          # Stage 4 lines: drug module + low-resource module
            line_rhs_multi!(du, x, xi, S.hb[1][j], L, P.p, (D=u, R=env.R), crowd, T)
        else
            w = L.input == :D ? T(u) : T(env.R)
            line_rhs!(du, x, xi, S.hb[1][j], L, P.p, w, crowd, T)
        end
    end
    nothing
end
function plant_f(x, u, P::Plant)
    T = promote_type(eltype(x), typeof(u))
    du = zeros(T, length(x)); plant_rhs!(du, x, u, P); du
end

# ── discrete transition F(x, u) over one decision interval (fixed-step RK4) ─────

const RK4_H = Float64(CFG["schedule"]["rk4_step_days"])
function F(x, u, P::Plant, Δ; h=RK4_H)
    n = max(1, ceil(Int, Δ / h - 1e-9)); h = Δ / n
    T = promote_type(eltype(x), typeof(u))
    y = T.(x)
    for _ in 1:n
        k1 = plant_f(y, u, P); k2 = plant_f(y .+ (h / 2) .* k1, u, P)
        k3 = plant_f(y .+ (h / 2) .* k2, u, P); k4 = plant_f(y .+ h .* k3, u, P)
        y = y .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    y
end
# dense trajectory of a schedule (u per interval) for plotting
function trajectory(x0, U, P::Plant, Δ; dt=0.1)
    ts = [0.0]; xs = [copy(x0)]; x = copy(x0); t = 0.0
    m = max(1, round(Int, Δ / dt))
    for u in U, _ in 1:m
        x = F(x, u, P, Δ / m); t += Δ / m
        push!(ts, t); push!(xs, x)
    end
    ts, xs
end

# Local linearisation around (x̄, ū) by forward-mode automatic differentiation:
#   A = ∂F/∂x, B = ∂F/∂u, one ForwardDiff Jacobian of F over z = [x; u] (nx + 1 dual directions).
# With ~9–12 states and one input, forward mode needs nx + 1 passes and reverse mode nx passes,
# so the cost is similar; forward mode differentiates through the RK4 loop directly and needs no tape.
# At a dose bound (u ≤ lo or u ≥ hi) B is the one-sided derivative into the feasible side (second-order
# one-sided difference): at u = 0 with no drug in the system the drug-effect states sit exactly on their
# max(E, 0) clamp, where AD takes the clamped branch and misses the dose sensitivity (2026-11-08).
function linearize(x, u, P::Plant, Δ; lo=-Inf, hi=Inf)
    nx = length(x)
    J = ForwardDiff.jacobian(z -> F(z[1:nx], z[nx+1], P, Δ), vcat(x, u))
    Fx = F(x, u, P, Δ); B = J[:, nx+1]
    s = u <= lo ? 1.0 : u >= hi ? -1.0 : 0.0
    s == 0 || (h = 1e-4 * s; B = (4 .* F(x, u + h, P, Δ) .- F(x, u + 2h, P, Δ) .- 3 .* Fx) ./ (2h))
    Fx, J[:, 1:nx], B
end

# Starting counts of every control run: the NOMINAL seeding of the coculture experiments
# (control_config.toml [initial]: 100 cells in total, split by the sensitive:resistant ratio, so
# "25:75" → N = 25, C = 75). Not a measured day-0 value: the zoo fits treat day 0 as the known seeding
# and the candidates' u0 as their starting state, so the controller starts from the design value too.
function nominal_seeding(ratio::AbstractString)
    s, r = parse.(Float64, split(ratio, ":"))
    tot = Float64(CFG["initial"]["seeding_total"])
    (N0=tot * s / (s + r), C0=tot * r / (s + r))
end

# initial state: the given counts, drug and toxicity states at 0
function x0_for(P::Plant, N0, C0)
    x = zeros(P.nx); x[P.iN] = N0; x[P.iC] = C0
    P.fp === nothing || (x[P.fp.iW] = 1.0)
    L = P.sys.lines
    L[2].action == :resource && error("resource states need initialisation")
    x
end
