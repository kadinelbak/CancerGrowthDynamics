# Two-input optimal control: cisplatin dose u (µM) and low resource / fasting v (0 = fed, 10% FBS;
# 1 = 0.5% FBS), with a host-weight cost and fasting protection of healthy tissue.
# Include after nothing else: it includes control_model.jl and control_mpc.jl.
#
# State x = [N, C, A_N,1..3, A_C,1..3, S_C,1..3, B, T, P, W]   (S_C: low-resource stress compartments of A2780cis)
# Control w = [u, v], piecewise constant over each 2-day decision interval.
#
# Cells: the scenario A/B normal-medium plant of control_model.jl (Stage 3 or monoculture drug response)
# with the fitted low-resource modules of each line added (Stage 1b/3b structure, Stage 4 strengths):
#   u_R = 0.955·v;  A2780: K_eff = K / (1 + c_N·ε_N(u)·u_R)
#                   A2780cis: extra loss (d_C^LR·ε_C(u)·u_R + d_C^acc·S_C,3)·C,  S_C chain driven by ε_C(u)·u_R
#   ε_X(u) = 1 + (ε_X^S4 − 1)·min(u/u_ref, 1)     (exact at both measured conditions; see control_config.toml)
# so v = 0 gives exactly the normal-medium plant, and (u = 1, v = 1) gives exactly the Stage 4 model.
# Host (control_config.toml [toxicity], [host]):
#   dB/dt = s_B·u − d_B·B
#   dT/dt = a_H·(1 − p_H·P)·(1 − e^(−B))·(1 − T) − r_H·T·(1 − T)      ([toxicity] model = "depillis")
#   dE/dt = k_N·Emax·γu/(EC50 + γu)·(1 − p_H·P) − k_N·E                ("xie2020", the default: kidney excess, Xie et al. 2020)
#   dP/dt = k_P·(v − P)
#   dW/dt = −k_W·v + r_W·(1 − W)·(1 − v)

include(joinpath(@__DIR__, "control_model.jl"))
include(joinpath(@__DIR__, "control_mpc.jl"))

const FCFG = CFG["fasting"]; const HCFG = CFG["host"]
default_host() = (k_W=Float64(HCFG["k_W"]), r_W=Float64(HCFG["r_W"]), k_P=Float64(HCFG["k_P"]), p_H=Float64(HCFG["p_H"]))

# The cell equations are simulated with the flattened, type-stable copy of control_model.jl (ModSpec,
# LineSpec, line_specs, spec_cells_rhs!), here with the low-resource modules and their free strengths ε.
# test_control_fasting.jl checks it against the zoo right-hand side (line_rhs_multi!) to 1e-12.

struct FPlant{TT<:NamedTuple,HT<:NamedTuple}
    name::String
    base::Plant                 # sys with the low-resource modules; p = [scenario params; Stage 4 params]
    iε::Vector{Int}             # positions in base.p of the low-resource strengths Stage 4 refits (none if it carries χ instead)
    ε::Vector{Float64}          # their Stage 4 (or bootstrap) values
    eps_mode::Symbol            # :interp or :stage4
    u_ref::Float64
    nx::Int; iN::Int; iC::Int; iB::Int; iT::Int; iP::Int; iW::Int
    K::Float64
    tox::TT; host::HT
    spec::Vector{LineSpec}
    koff::Float64               # fade rate of low-resource stress after refeeding; NaN = the fitted build-up rate (the zoo model)
end

shift_slots(s, np) = NamedTuple{keys(s)}(Tuple(v.i > 0 ? Slot(v.i + np, v.v) : v for v in values(s)))

function fast_plant(scen::Symbol; b::Union{Nothing,Int}=nothing, tox=default_tox(), host=default_host(),
                    eps_mode=Symbol(FCFG["eps_mode"]), koff=NaN)
    P0 = scenario_plant(scen; b=b, tox=tox)
    P0 === nothing && return nothing
    # Stage 4 adds each line's low-resource module and whatever it refits: low-resource strengths ε (free ε slots
    # of the modules) and/or a drug × low-resource interaction χ (the line's chi slot); its parameters are
    # appended to the scenario's
    c4 = zgroup("s4").cands[1]
    hasproperty(c4, :slots) && :K in c4.slots && error("a Stage 4 model that refits K is not supported by the control plant")
    p4 = b === nothing ? copy(c4.params) : ZOO_STATE.results["s4"].res.pars[1][b, :]
    any(isnan, p4) && return nothing
    np = length(P0.p)
    shift(sl) = sl.i > 0 ? Slot(sl.i + np, sl.v) : sl
    lines = map(1:2) do j
        L4 = c4.sys.lines[j]
        L = merge(P0.sys.lines[j], (extra=[merge(M, (s=shift_slots(M.s, np),)) for M in L4.extra],))
        hasproperty(L4, :chi) ? merge(L, (chi=shift(L4.chi),)) : L
    end
    p = vcat(P0.p, p4)
    iε = Int[np + M.s.ε.i for j in 1:2 for M in c4.sys.lines[j].extra if M.s.ε.i > 0]
    ε = p[iε]
    base = make_plant(P0.name, lines, P0.sys.alpha, p; uR=0.0, tox=tox, fast=false)     # the FPlant has its own spec
    n = base.sys.n
    name = (scen == :A ? "A: Stage 3 coculture model as selected" : "B: monoculture drug response for both lines") *
           " + fitted low-resource modules"
    FPlant(name, base, iε, ε, eps_mode, Float64(FCFG["eps_u_ref"]),
           n + 4, base.iN, base.iC, n + 1, n + 2, n + 3, n + 4, base.K, tox, host, line_specs(base, iε), Float64(koff))
end

# The 14-day controller's low-resource environment (control_model.jl Plant interface): the same plant
# held at v = 1 for the whole run, with no fasting protection of healthy tissue (p_H = 0), because low
# resource is the culture's environment there, not a fast the host chooses. ε(u) follows eps_mode as above,
# so u = 0 gives the low-resource-only model and u = 1 µM the Stage 4 model.
function lr_plant(scen::Symbol; b::Union{Nothing,Int}=nothing, tox=default_tox())
    FP = fast_plant(scen; b=b, tox=tox, host=merge(default_host(), (p_H=0.0,)))
    FP === nothing && return nothing
    B = FP.base
    name = (scen == :A ? "A (low resource): Stage 3 coculture model" : "B (low resource): monoculture drug response") *
           " + fitted low-resource modules, dose-dependent strength ε(u)"
    Plant(name, B.sys, B.p, FP.nx, FP.iN, FP.iC, FP.iB, FP.iT, FP.K, tox, U_R, FP, nothing)
end

εfactor(FP::FPlant, u) = FP.eps_mode == :stage4 ? one(u) : min(u / FP.u_ref, one(u))

function fast_rhs_ref!(du, x, w, FP::FPlant)       # reference: the zoo right-hand side
    u, v = w[1], w[2]
    T = promote_type(eltype(x), eltype(w))
    B = FP.base; S = B.sys; env = S.envs[1]
    p = T.(B.p); f = εfactor(FP, T(u))
    for (j, i) in enumerate(FP.iε)
        p[i] = 1 + (FP.ε[j] - 1) * f
    end
    envc = (D=T(u), R=T(U_R) * v)
    for (j, l) in enumerate(env.lines)
        xi = S.xi[1][j]; L = S.lines[l]
        crowd = max(x[xi+1], zero(T)); other = zero(T)
        for (m, lm) in enumerate(env.lines)
            m == j && continue
            y = max(x[S.xi[1][m]+1], zero(T))
            crowd += sv(S.alpha[l][lm], p, T) * y; other += y
        end
        line_rhs_multi!(du, x, xi, S.hb[1][j], L, p, envc, crowd, T; other=other)
    end
    Bd = x[FP.iB]; Tx = x[FP.iT]; Pr = x[FP.iP]; Wt = x[FP.iW]; tp = FP.tox; h = FP.host
    du[FP.iB] = tp.s_B * u - tp.d_B * Bd
    du[FP.iT] = tp.xie ? xie_drive(tp, u) * (1 - h.p_H * Pr) - tp.kN * Tx :
                         tp.a_H * (1 - h.p_H * Pr) * (1 - exp(-max(Bd, zero(T)))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
    du[FP.iP] = h.k_P * (v - Pr)
    du[FP.iW] = -h.k_W * v + h.r_W * (1 - Wt) * (1 - v)
    nothing
end
function fast_rhs!(du, x, w, FP::FPlant)
    u, v = w[1], w[2]
    T = promote_type(eltype(x), eltype(w)); z = zero(T)
    spec_cells_rhs!(du, x, T(u), T(U_R) * v, εfactor(FP, T(u)), FP.spec, FP.ε, FP.koff, v)     # control_model.jl
    Bd = x[FP.iB]; Tx = x[FP.iT]; Pr = x[FP.iP]; Wt = x[FP.iW]; tp = FP.tox; h = FP.host
    du[FP.iB] = tp.s_B * u - tp.d_B * Bd
    du[FP.iT] = tp.xie ? xie_drive(tp, u) * (1 - h.p_H * Pr) - tp.kN * Tx :
                         tp.a_H * (1 - h.p_H * Pr) * (1 - exp(-max(Bd, z))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
    du[FP.iP] = h.k_P * (v - Pr)
    du[FP.iW] = -h.k_W * v + h.r_W * (1 - Wt) * (1 - v)
    nothing
end
fast_f(x, w, FP::FPlant) = (T = promote_type(eltype(x), eltype(w)); du = zeros(T, length(x)); fast_rhs!(du, x, w, FP); du)

function F2(x, w, FP::FPlant, Δ; h=RK4_H)
    n = max(1, ceil(Int, Δ / h - 1e-9)); h = Δ / n
    y = promote_type(eltype(x), eltype(w)).(x)
    for _ in 1:n
        k1 = fast_f(y, w, FP); k2 = fast_f(y .+ (h / 2) .* k1, w, FP)
        k3 = fast_f(y .+ (h / 2) .* k2, w, FP); k4 = fast_f(y .+ h .* k3, w, FP)
        y = y .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    y
end
function trajectory2(x0, U, FP::FPlant, Δ; dt=0.1)
    ts = [0.0]; xs = [copy(x0)]; x = copy(x0); t = 0.0; m = max(1, round(Int, Δ / dt))
    for k in axes(U, 1), _ in 1:m
        x = F2(x, U[k, :], FP, Δ / m); t += Δ / m
        push!(ts, t); push!(xs, x)
    end
    ts, xs
end
# A_k = ∂F/∂x (nx×nx), B_k = ∂F/∂w (nx×2), one forward-mode Jacobian over [x; u; v]
# at an input bound (w ≤ lo or w ≥ hi) that column of B is the one-sided derivative into the feasible side,
# as in linearize (control_model.jl): AD misses the dose sensitivity on the max(E, 0) clamp
function linearize2(x, w, FP::FPlant, Δ; lo=fill(-Inf, length(w)), hi=fill(Inf, length(w)))
    nx = length(x)
    J = ForwardDiff.jacobian(z -> F2(z[1:nx], z[nx+1:end], FP, Δ), vcat(x, w))
    Fx = F2(x, w, FP, Δ); B = J[:, nx+1:end]
    for i in eachindex(w)
        s = w[i] <= lo[i] ? 1.0 : w[i] >= hi[i] ? -1.0 : 0.0
        s == 0 && continue
        e = (1e-4 * s) .* (eachindex(w) .== i)
        B[:, i] = (4 .* F2(x, w .+ e, FP, Δ) .- F2(x, w .+ 2 .* e, FP, Δ) .- 3 .* Fx) ./ (2e-4 * s)
    end
    Fx, J[:, 1:nx], B
end
function x0_fast(FP::FPlant, N0, C0)
    x = zeros(FP.nx); x[FP.iN] = N0; x[FP.iC] = C0; x[FP.iW] = 1.0
    x
end

# ── cost ─────────────────────────────────────────────────────────────────────
struct Weights2
    qB::Float64; qR::Float64; qT::Float64; qW::Float64; ru::Float64; rv::Float64; tf::Float64
    μC::Float64; μT::Float64; μW::Float64; δ::Float64
    Bc::Float64; Btox::Float64; Wmin::Float64; umax::Float64; K::Float64
    qP::Float64; Bp::Float64; sP::Float64; qRT::Float64      # time to progression (control_mpc.jl)
    BpN::Float64; BpC::Float64                                # per-line thresholds of the "lines" rule (control_mpc.jl line_thresholds)
    sB::Float64                                               # scale of the inside of the host barriers (control_mpc.jl host_barrier)
end
# untreated, unfasted sensitive-only plateau of the fasting plant (normal medium)
function sensitive_plateau(FP::FPlant; days=400.0)
    x = x0_fast(FP, 0.05FP.K, 0.0)
    for _ in 1:round(Int, days / 2); x = F2(x, [0.0, 0.0], FP, 2.0); end
    x[FP.iN]
end
# bc_frac: B_cancer as a fraction of K (0 = no soft burden penalty); tf: terminal factor (both swept)
function weights2(FP::FPlant; umax=Float64(CFG["dose"]["u_max_data"]), Wmin=Float64(HCFG["W_min"]), w=nothing, rv=Float64(FCFG["rv"]),
                  bc_frac=Float64(CFG["cancer"]["B_cancer_fraction_of_K"]), tf=Float64(CFG["weights"]["terminal_factor"]),
                  prog_frac=Float64(CFG["progression"]["threshold_fraction_of_K"]), ts=tiebreak_scale())
    c = CFG["weights"]
    qB, qR, qT, ru, qP, qRT = objective_weights(w; ts=ts)
    Bp, sP = progression_scale(FP.K; prog_frac=prog_frac)
    BpN, BpC = prog_rule() === :lines && prog_frac > 0 ? line_thresholds(FP.K, sensitive_plateau(FP); prog_frac=prog_frac) : (Inf, Inf)
    sc = qP > 0 ? ts : 1.0           # host weight and fasting costs are tie-breakers too (control_mpc.jl)
    Weights2(qB, qR, qT, HCFG["qW"] * sc, ru, rv * sc, tf, c["mu_cancer"], c["mu_tox"], HCFG["mu_host"], c["barrier_delta"],
             bc_frac > 0 ? bc_frac * FP.K : Inf, tox_limit(FP.tox), Wmin, umax, FP.K, qP, Bp, sP, qRT, BpN, BpC, sc)
end
wloss(x, FP, W) = (1 - x[FP.iW]) / (1 - W.Wmin)
# the soft burden penalty at B_cancer follows control_config [cancer] barrier: "none" (default since
# 2026-11-01, as B_cancer = B_prog and the time-to-progression term already prices it) leaves it out
const BURDEN_BARRIER2 = String(get(CFG["cancer"], "barrier", "relaxed_log")) != "none"
function stage_cost2(x, w, FP::FPlant, W::Weights2)
    N, C, T = x[FP.iN], x[FP.iC], x[FP.iT]
    progression_cost(N, C, W) + burden_cost(N, C, W) + W.qT * (T / W.Btox)^2 + W.qW * wloss(x, FP, W)^2 +
        W.ru * (w[1] / W.umax)^2 + W.rv * w[2]^2 +
        (BURDEN_BARRIER2 ? relaxed_barrier(1 - (N + C) / W.Bc, W.μC, W.δ) : zero(N)) + host_barrier(1 - T / W.Btox, W.μT, W.δ, W.sB) +
        host_barrier(1 - wloss(x, FP, W), W.μW, W.δ, W.sB)
end
terminal_cost2(x, FP::FPlant, W::Weights2) = W.tf * (progression_cost(x[FP.iN], x[FP.iC], W) + burden_cost(x[FP.iN], x[FP.iC], W) +
                                                     W.qT * (x[FP.iT] / W.Btox)^2 + W.qW * wloss(x, FP, W)^2) +
    W.qRT * max(x[FP.iC], zero(x[FP.iC])) / W.K

# blocked (dose, fast) plan → one row per interval: the dose is held over the block, a fast only in its first interval
expand_blocks2(Ub, b) = reduce(vcat, [[Ub[j, 1] (i == 1 ? Ub[j, 2] : 0.0)] for j in eachindex(b) for i in 1:b[j]])
function total_cost2(x0, U, FP, Δ, W)
    x = copy(x0); J = 0.0
    for k in axes(U, 1)
        x = F2(x, U[k, :], FP, Δ); J += stage_cost2(x, U[k, :], FP, W)
    end
    J + terminal_cost2(x, FP, W)
end
# the cost the exact search minimises: as total_cost2, with the progression term averaged over PROG_SUB pieces
function plan_cost2(x0, U, FP, Δ, W)
    x = copy(x0); J = 0.0
    for k in axes(U, 1)
        w = U[k, :]; pc = 0.0
        for _ in 1:PROG_SUB
            x = F2(x, w, FP, Δ / PROG_SUB); pc += progressed(x[FP.iN], x[FP.iC], W)
        end
        J += stage_cost2(x, w, FP, W) - progression_cost(x[FP.iN], x[FP.iC], W) + W.qP * pc / PROG_SUB
    end
    J + terminal_cost2(x, FP, W)
end

# The stage cost split into its terms, summed over the intervals of an applied schedule (no terminal
# term), as cost_terms in control_mpc.jl plus the host terms of fasting: total = Σ_k stage_cost2.
function cost_terms2(FP::FPlant, x0, U, Δ, W::Weights2)
    x = copy(x0); z = zeros(10)
    for k in axes(U, 1)
        w = U[k, :]; x = F2(x, w, FP, Δ); N, C, T = x[FP.iN], x[FP.iC], x[FP.iT]
        z .+= [W.qB * (max(N, 0.0) + max(C, 0.0)) / W.K, W.qR * max(C, 0.0) / W.K, W.qT * (T / W.Btox)^2, W.ru * (w[1] / W.umax)^2,
               BURDEN_BARRIER2 ? relaxed_barrier(1 - (N + C) / W.Bc, W.μC, W.δ) : 0.0, host_barrier(1 - T / W.Btox, W.μT, W.δ, W.sB),
               W.qW * wloss(x, FP, W)^2, W.rv * w[2]^2, host_barrier(1 - wloss(x, FP, W), W.μW, W.δ, W.sB), progression_cost(N, C, W)]
    end
    (burden=z[1], resistance=z[2], toxicity=z[3], dose=z[4], penalty_burden=z[5], barrier_tox=z[6],
     weight=z[7], fasting=z[8], barrier_weight=z[9], progression=z[10], total=sum(z))
end

# adjoint (vector–Jacobian) gradient and Gauss–Newton Hessian of J over the stacked controls
function grad_hess2(x0, U, FP::FPlant, Δ, W::Weights2)
    H = size(U, 1); m = size(U, 2)
    xs = [copy(x0)]; As = Matrix{Float64}[]; Bs = Matrix{Float64}[]
    for k in 1:H
        xn, A, B = linearize2(xs[k], U[k, :], FP, Δ; lo=[0.0, 0.0], hi=[W.umax, 1.0])
        push!(xs, xn); push!(As, A); push!(Bs, B)
    end
    # adjoint gradient
    g = zeros(H * m)
    μ = ForwardDiff.gradient(x -> stage_cost2(x, U[H, :], FP, W), xs[H+1]) .+ ForwardDiff.gradient(x -> terminal_cost2(x, FP, W), xs[H+1])
    for k in H:-1:1
        gw = ForwardDiff.gradient(ww -> stage_cost2(xs[k+1], ww, FP, W), U[k, :])
        g[(k-1)*m+1:k*m] = gw .+ Bs[k]' * μ
        k > 1 && (μ = ForwardDiff.gradient(x -> stage_cost2(x, U[k-1, :], FP, W), xs[k]) .+ As[k]' * μ)
    end
    # Gauss–Newton Hessian
    nx = length(x0); S = zeros(nx, H * m); Hs = zeros(H * m, H * m)
    for k in 1:H
        S = As[k] * S; S[:, (k-1)*m+1:k*m] .+= Bs[k]
        Q = ForwardDiff.hessian(x -> stage_cost2(x, U[k, :], FP, W), xs[k+1])
        k == H && (Q = Q + ForwardDiff.hessian(x -> terminal_cost2(x, FP, W), xs[k+1]))
        Hs .+= S' * psd(Q) * S
        Hs[(k-1)*m+1:k*m, (k-1)*m+1:k*m] .+= psd(ForwardDiff.hessian(ww -> stage_cost2(xs[k+1], ww, FP, W), U[k, :]))
    end
    g, (Hs + Hs') / 2 + 1e-9 * LA.I
end

# ── SQP (continuous u ∈ [0, u_max], v ∈ [0, 1]); variables ordered (k, input) → (k−1)·2 + i ─────
function sqp_solve2(x0, U0, FP::FPlant, Δ, W::Weights2; max_iter=Int(CFG["sqp"]["max_iter"]), tol=Float64(CFG["sqp"]["tol"]))
    H = size(U0, 1); m = 2
    lo = zeros(H * m); hi = repeat([W.umax, 1.0], H)
    z = clamp.(vec(permutedims(U0)), lo, hi)
    Um(z) = permutedims(reshape(z, m, H))
    J = total_cost2(x0, Um(z), FP, Δ, W); it = 0; Jold = J
    for outer it in 1:max_iter
        g, Hs = grad_hess2(x0, Um(z), FP, Δ, W)
        d = box_qp(Hs, g .- Hs * z, lo, hi, z) .- z
        maximum(abs.(d)) < 1e-8 && break
        α = 1.0; slope = LA.dot(g, d); accepted = false
        for _ in 1:25
            zn = z .+ α .* d; Jn = total_cost2(x0, Um(zn), FP, Δ, W)
            if Jn <= J + 1e-4 * α * min(slope, 0.0)
                z, Jold, J = zn, J, Jn; accepted = true; break
            end
            α /= 2
        end
        if !accepted         # QP step rejected: a projected-gradient step instead (control_mpc.jl pg_step)
            pg = pg_step(zz -> total_cost2(x0, Um(zz), FP, Δ, W), J, g, z, lo, hi, max(LA.opnorm(Hs), 1e-12))
            pg === nothing || ((z, Jold, J) = (pg[1], J, pg[2]); accepted = true)
        end
        (!accepted || abs(Jold - J) < tol * (1 + abs(J))) && break
    end
    (U=Um(z), J=J, iters=it)
end

# allocation-free RK4 transition for the exact search (Float64 only; same steps as F2)
struct RK4Buf
    k1::Vector{Float64}; k2::Vector{Float64}; k3::Vector{Float64}; k4::Vector{Float64}; tmp::Vector{Float64}
end
RK4Buf(n) = RK4Buf(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))
function F2!(y, x, w, FP::FPlant, Δ, b::RK4Buf; h=RK4_H)
    n = max(1, ceil(Int, Δ / h - 1e-9)); h = Δ / n
    copyto!(y, x)
    for _ in 1:n
        fast_rhs!(b.k1, y, w, FP); @. b.tmp = y + (h / 2) * b.k1
        fast_rhs!(b.k2, b.tmp, w, FP); @. b.tmp = y + (h / 2) * b.k2
        fast_rhs!(b.k3, b.tmp, w, FP); @. b.tmp = y + h * b.k3
        fast_rhs!(b.k4, b.tmp, w, FP); @. y += (h / 6) * (b.k1 + 2 * b.k2 + 2 * b.k3 + b.k4)
    end
    y
end

# ── exact search over measured doses × {fed, fasting} with the fasting-schedule rules ───────────────
# Depth-first over every admissible sequence of (dose, fed/fasting) pairs. Exact: a branch is only cut
# when its partial cost already exceeds the best complete schedule (every stage cost is ≥ 0). To cut
# early, the search starts from the cost of the warm-start plan (the previous plan shifted by one
# interval) and expands the children of each node cheapest-first.
# run = consecutive fasting intervals ending at k−1; gap = fed intervals since the last fast.
admissible(v, run, gap, maxrun, mingap) = v == 0 || (run < maxrun && (run > 0 || gap >= mingap))
# H: an Int (single intervals) or the look-ahead blocks (control_mpc.jl). In a block the dose holds for
# every interval and a fast applies to the block's FIRST interval only (fed afterwards), so the fasting
# rules (48-h fasts, refeeding between) still hold; the cost is summed over every interval.
function discrete_solve2(x0, FP::FPlant, Δ, H, W::Weights2; ulevels=Float64.(CFG["dose"]["levels_data"]),
                         vlevels=Float64.(FCFG["levels"]), maxrun=Int(FCFG["max_consecutive_intervals"]),
                         mingap=Int(FCFG["min_gap_intervals"]), run0=0, gap0=10^6, Uinit=nothing)
    blk = as_blocks(H); H = length(blk)
    acts = [(u, v) for v in vlevels for u in ulevels]; na = length(acts); nx = length(x0)
    buf = RK4Buf(nx); tmp = zeros(nx)
    # one block from x into y: returns the block's cost and the fasting run/gap after it
    function block!(y, x, u, v, k, run, gap)
        local cb = 0.0; copyto!(tmp, x)          # `local`: not the warm-start block's `c` (closure capture)
        for j in 1:blk[k]
            vj = j == 1 ? v : 0.0
            # time spent progressed, averaged over PROG_SUB pieces of the interval (control_mpc.jl: step_progressed)
            local pc = 0.0
            for _ in 1:PROG_SUB
                F2!(y, tmp, (u, vj), FP, Δ / PROG_SUB, buf); copyto!(tmp, y); pc += progressed(y[FP.iN], y[FP.iC], W)
            end
            cb += stage_cost2(y, (u, vj), FP, W) - progression_cost(y[FP.iN], y[FP.iC], W) + W.qP * pc / PROG_SUB
            run, gap = vj > 0 ? (run + 1, 0) : (0, gap + 1)
        end
        cb, run, gap
    end
    X = [[zeros(nx) for _ in 1:na] for _ in 1:H]; cst = [zeros(na) for _ in 1:H]; ord = [collect(1:na) for _ in 1:H]
    best = Ref(Inf); bestU = zeros(H, 2); U = zeros(H, 2)
    if Uinit !== nothing                      # upper bound from the warm-start plan, if it obeys the rules
        x = copy(x0); c = 0.0; run, gap = run0, gap0; ok = true; y = zeros(nx)
        for k in 1:H
            u, v = Uinit[k, 1], Uinit[k, 2]
            admissible(v, run, gap, maxrun, mingap) || (ok = false; break)
            ck, run, gap = block!(y, x, u, v, k, run, gap); x, y = y, x; c += ck
        end
        ok && (best[] = c + terminal_cost2(x, FP, W); bestU .= Uinit)
    end
    RG = [zeros(Int, na, 2) for _ in 1:H]            # fasting run/gap after each child block
    function dfs(x, k, cost, run, gap)
        if k > H
            local cl = cost + terminal_cost2(x, FP, W)
            cl < best[] && (best[] = cl; bestU .= U)
            return
        end
        for a in 1:na
            u, v = acts[a]
            if admissible(v, run, gap, maxrun, mingap)
                ck, r2, g2 = block!(X[k][a], x, u, v, k, run, gap); cst[k][a] = cost + ck
                RG[k][a, 1] = r2; RG[k][a, 2] = min(g2, 10^6)
            else
                cst[k][a] = Inf
            end
        end
        sortperm!(ord[k], cst[k])
        for a in ord[k]
            cst[k][a] >= best[] && break
            u, v = acts[a]; U[k, 1] = u; U[k, 2] = v
            dfs(X[k][a], k + 1, cst[k][a], RG[k][a, 1], RG[k][a, 2])
        end
    end
    dfs(copy(x0), 1, 0.0, run0, gap0)
    (U=copy(bestU), J=best[])
end

# Receding horizon on plant Pp, planning with Pm; the fasting rules carry over between decisions.
# Fixed look-ahead H at every decision, past the end of the run if need be (shrink=true: comparison only),
# so the schedule is the first decisions of a continuing policy. The planner observes only N and C
# (× noise[k, :] when given; control_mpc.jl: observe!) and reconstructs the drug-effect and stress chains,
# B, T, P and W by simulating its own model under the applied doses and fasting.
# H: an Int (single intervals) or the look-ahead blocks (discrete mode); continuous mode plans every interval.
# Early stop as closed_loop (control_mpc.jl: stop_now): Bp, confirm, settle_day.
function closed_loop2(Pm::FPlant, Pp::FPlant, x0, Δ, n_dec, H, W::Weights2; mode=:discrete, shrink=false, noise=nothing,
                      Bp=Inf, confirm=0, settle_day=Inf, kw...)
    blk = mode == :continuous ? fill(1, sum(as_blocks(H))) : as_blocks(H); nb = length(blk)
    x = copy(x0); xm = copy(x0); applied = zeros(n_dec, 2); run = 0; gap = 10^6
    Uwarm = repeat([mode == :discrete ? 0.0 : 0.5  0.0], nb)
    above = 0; stop_day = NaN; stop_reason = ""
    for k in 1:n_dec
        observe!(xm, x, Pm, Pp, noise, k)
        bk = shrink ? cut_blocks(blk, n_dec - k + 1) : blk; Hk = length(bk)
        if mode == :discrete
            sol = discrete_solve2(xm, Pm, Δ, bk, W; run0=run, gap0=gap, Uinit=Uwarm[1:Hk, :], kw...)
        else
            # SQP is local and the progression step makes the cost non-convex: start from the cheaper of the
            # shifted warm start and the exact measured-level plan (blocked as in discrete mode), so the
            # continuous plan is never worse than it (SQP does not increase the cost)
            bd = shrink ? cut_blocks(as_blocks(H), n_dec - k + 1) : as_blocks(H)
            U0 = Uwarm[1:Hk, :]; Ud = expand_blocks2(discrete_solve2(xm, Pm, Δ, bd, W; run0=run, gap0=gap).U, bd)
            total_cost2(xm, Ud, Pm, Δ, W) < total_cost2(xm, U0, Pm, Δ, W) && (U0 = Ud)
            sol = sqp_solve2(xm, U0, Pm, Δ, W)
        end
        w = sol.U[1, :]; applied[k, :] = w
        run = w[2] > 0 ? run + 1 : 0; gap = w[2] > 0 ? 0 : gap + 1
        x = F2(x, w, Pp, Δ); xm = F2(xm, w, Pm, Δ)
        Uwarm = vcat(sol.U[2:end, :], repeat([sol.U[end, 1] 0.0], nb - Hk + 1))[1:nb, :]     # shift; pad fed (a bound only)
        z = isfinite(Bp) ? prog_index(x[Pp.iN], x[Pp.iC], W) : 0.0; above = z >= 1 ? above + 1 : 0
        r = k < n_dec ? stop_now(above, z, k * Δ, 1.0, confirm, settle_day) : nothing
        if r !== nothing
            # the rest follows the last plan, then its last dose, fed (no plan is computed after the stop)
            plan = mode == :discrete ? expand_blocks2(sol.U, bk) : sol.U
            for j in 1:(n_dec - k)
                wj = j + 1 <= size(plan, 1) ? plan[j+1, :] : [plan[end, 1], 0.0]
                applied[k+j, :] = wj; x = F2(x, wj, Pp, Δ)
            end
            stop_day = k * Δ; stop_reason = r
            break
        end
    end
    (U=applied, x=x, stop_day=stop_day, stop_reason=stop_reason)
end

# ── baseline schedules ────────────────────────────────────────────────────────
function run_policy2(policy, FP::FPlant, x0, Δ, n)
    x = copy(x0); U = zeros(n, 2)
    for k in 1:n
        U[k, :] = policy(x, k); x = F2(x, U[k, :], FP, Δ)
    end
    U
end
function threshold_policy2(FP::FPlant, W::Weights2; dose=1.0, hi=0.5, lo=0.25)
    on = Ref(false)
    (x, k) -> begin
        tot = x[FP.iN] + x[FP.iC]
        tot >= hi * W.Bc && (on[] = true); tot <= lo * W.Bc && (on[] = false)
        [on[] ? dose : 0.0, 0.0]
    end
end
baseline_policies2(FP::FPlant, W::Weights2) = [
    ("No treatment", (x, k) -> [0.0, 0.0]),
    ("Constant 1 µM", (x, k) -> [1.0, 0.0]),
    ("Maximum measured dose (1.47 µM)", (x, k) -> [1.47, 0.0]),
    ("Periodic 1 µM, 2 days on / 2 off", (x, k) -> [isodd(k) ? 1.0 : 0.0, 0.0]),
    ("Threshold adaptive (Zhang-style)", threshold_policy2(FP, W)),
    ("Fasting only, 48 h every 6 days", (x, k) -> [0.0, k % 3 == 1 ? 1.0 : 0.0]),
    ("Periodic 1 µM with a 48-h fast on each dosing interval", (x, k) -> isodd(k) ? [1.0, 1.0] : [0.0, 0.0]),
    ("Constant 1 µM + 48-h fast every 6 days", (x, k) -> [1.0, k % 3 == 1 ? 1.0 : 0.0]),
    # a fixed fasting schedule that keeps W ≥ W_min (user, 2026-11-12: the two above break the weight limit)
    ("Constant 1 µM + 48-h fast every 8 days", (x, k) -> [1.0, k % 4 == 1 ? 1.0 : 0.0]),
]

function outcomes2(FP::FPlant, x0, U, Δ, W::Weights2; t_stop=Inf)
    ts, xs = trajectory2(x0, U, FP, Δ)
    N = [x[FP.iN] for x in xs]; C = [x[FP.iC] for x in xs]; T = [x[FP.iT] for x in xs]
    Pr = [x[FP.iP] for x in xs]; Wt = [x[FP.iW] for x in xs]; tot = N .+ C
    over = findfirst(>=(W.Bc), tot)
    days_over = sum((ts[i+1] - ts[i] for i in 1:length(ts)-1 if tot[i+1] >= W.Bc); init=0.0)     # days with N + C ≥ B_cancer
    days_half_K = sum((ts[i+1] - ts[i] for i in 1:length(ts)-1 if tot[i+1] >= W.K / 2); init=0.0)
    tp, ip = progression_window(ts, prog_index.(N, C, Ref(W)), 1.0, t_stop)       # up to the closed loop's stop (control_mpc.jl)
    (N_final=N[end], C_final=C[end], total_final=tot[end], resistant_fraction=C[end] / max(tot[end], 1e-9),
     ttp_first=tp.first, ttp=tp.loss, days_progressed=tp.days, N_prog=ip > 0 ? N[ip] : NaN, C_prog=ip > 0 ? C[ip] : NaN,
     stop_day=isfinite(t_stop) ? t_stop : NaN,
     max_total=maximum(tot), max_tox=maximum(T), min_weight=minimum(Wt), cum_dose=sum(U[:, 1]) * Δ, fast_days=sum(U[:, 2]) * Δ,
     cancer_violation=over !== nothing, tox_violation=maximum(T) >= W.Btox, weight_violation=minimum(Wt) <= W.Wmin, days_over=days_over, days_half_K=days_half_K,
     ts=ts, N=N, C=C, T=T, P=Pr, W=Wt)
end
