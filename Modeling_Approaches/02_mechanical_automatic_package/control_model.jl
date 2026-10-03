# Control plant for optimal control of cisplatin dosing, built from the selected compartment models.
#
# The biology is NOT refitted here. The cell equations are the zoo's own right-hand side
# (zoo_models.jl: line_rhs!) with the parameters of the carried-forward models in
# outputs/zoo_fits/zoo_state.jls. Only the dose input and a toxicity proxy are added:
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
#       Stage 2 or refitted no lower than half of it, never switched off: user decision 2026-10-01)
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

# resistant kill β_C of the carried Stage 3 model, and where it comes from
function stage3_betaC()
    c3 = zgroup("s3").cands[1]; s = c3.sys.lines[2].s.β
    s.i > 0 ? (value=c3.params[s.i], free=true, source="Stage 3 (refitted, floor $(round(c3.bounds[s.i][1]; sigdigits=3)) /day)") :
              (value=s.v, free=false, source="Stage 2 (inherited)")
end

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
end

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

function make_plant(name, lines, alpha, p; uR=0.0, tox=default_tox())
    sys = build_system(lines, [(D=0.0, R=uR, lines=[1, 2])], alpha)
    n = sys.n
    Kslot = lines[1].s.K
    K = Kslot.i > 0 ? p[Kslot.i] : Kslot.v
    Plant(name, sys, p, n + 2, sys.xi[1][1] + 1, sys.xi[1][2] + 1, n + 1, n + 2, K, tox, uR, nothing)
end

# The Stage 3 coculture lines and parameters for a scenario (point estimates or a bootstrap draw b)
function scenario_plant(scen::Symbol; b::Union{Nothing,Int}=nothing, tox=default_tox())
    g3 = zgroup("s3"); c3 = g3.cands[1]
    LN, LC = c3.sys.lines
    draw(gid, i) = (P = ZOO_STATE.results[gid].res.pars[1]; b === nothing ? zgroup(gid).cands[1].params[i] : P[b, i])
    fixed_N = Dict{Symbol,Float64}(); fixed_C = Dict{Symbol,Float64}()
    if b !== nothing       # carried parameters from the same bootstrap replicate of their own stage
        fixed_N[:r] = draw("s1-a2780", 1); fixed_C[:r] = draw("s1-a2780cis", 1); fixed_C[:θ] = draw("s1-a2780cis", 3)
        fixed_N[:h] = draw("s2-a2780", 2); fixed_N[:k] = draw("s2-a2780", 3)
        fixed_C[:h] = draw("s2-a2780cis", 2); fixed_C[:k] = draw("s2-a2780cis", 3)
    end
    if scen == :A
        p = b === nothing ? copy(c3.params) : ZOO_STATE.results["s3"].res.pars[1][b, :]
        # an inherited resistant kill follows its own Stage 2 bootstrap draw, as h and k do
        b !== nothing && LC.s.β.i == 0 && (fixed_C[:β] = draw("s2-a2780cis", 1))
        name = "A: Stage 3 coculture model as selected"
    else
        ce0 = AN_STATE[:validation].ce1[1].params          # K, α_NC, α_CN, d_N fitted to untreated cocultures
        p = b === nothing ? copy(ce0) : ZOO_STATE.results["s3"].res.pars[1][b, 1:4] .* (ce0 ./ c3.params[1:4])
        fixed_N[:β] = draw("s2-a2780", 1); fixed_C[:β] = draw("s2-a2780cis", 1)
        name = "B: monoculture drug response for both lines, crowding from untreated cocultures"
    end
    any(isnan, p) && return nothing
    make_plant(name, [set_slots(LN, fixed_N), set_slots(LC, fixed_C)], c3.sys.alpha, p; tox=tox)
end

# Low-resource environment: lr_plant(scen) is defined in control_fasting.jl (the two-input plant held at
# v = 1), so the 14-day and 30-day controllers use the same dose-dependent drug × low-resource model.

# continuous-time right-hand side for a given dose u (µM)
function plant_rhs!(du, x, u, P::Plant)
    T = promote_type(eltype(x), typeof(u))
    P.fp === nothing || return fast_rhs!(du, x, (T(u), one(T)), P.fp)      # low resource: fasting plant at v = 1
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
    B = x[P.iB]; Tx = x[P.iT]; tp = P.tox
    du[P.iB] = tp.s_B * u - tp.d_B * B
    du[P.iT] = tp.a_H * (1 - exp(-max(B, zero(T)))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
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
function linearize(x, u, P::Plant, Δ)
    nx = length(x)
    J = ForwardDiff.jacobian(z -> F(z[1:nx], z[nx+1], P, Δ), vcat(x, u))
    F(x, u, P, Δ), J[:, 1:nx], J[:, nx+1]
end

# initial state: observed counts, drug and toxicity states at 0
function x0_for(P::Plant, N0, C0)
    x = zeros(P.nx); x[P.iN] = N0; x[P.iC] = C0
    P.fp === nothing || (x[P.fp.iW] = 1.0)
    L = P.sys.lines
    L[2].action == :resource && error("resource states need initialisation")
    x
end
