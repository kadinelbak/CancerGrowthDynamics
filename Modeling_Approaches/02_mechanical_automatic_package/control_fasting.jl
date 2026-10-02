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
#   dT/dt = a_H·(1 − p_H·P)·(1 − e^(−B))·(1 − T) − r_H·T·(1 − T)
#   dP/dt = k_P·(v − P)
#   dW/dt = −k_W·v + r_W·(1 − W)·(1 − v)

include(joinpath(@__DIR__, "control_model.jl"))
include(joinpath(@__DIR__, "control_mpc.jl"))

const FCFG = CFG["fasting"]; const HCFG = CFG["host"]
default_host() = (k_W=Float64(HCFG["k_W"]), r_W=Float64(HCFG["r_W"]), k_P=Float64(HCFG["k_P"]), p_H=Float64(HCFG["p_H"]))

# Flattened, type-stable copy of the plant's cell equations (the same modules and values as the zoo's
# line_rhs_multi!, which is generic and dispatches dynamically, about 50× slower). test_control.jl
# checks it against the zoo right-hand side to 1e-12. Only the module types of the carried models
# are supported; anything else stops with an error rather than silently differing.
struct ModSpec
    input_D::Bool; kernel_chain::Int; k::Float64; o::Int     # chain length (0 = immediate), rate, first hidden index − 1
    action::Symbol; β::Float64; h::Float64; c::Float64; d0::Float64; ec50::Float64
    ε::Float64; εj::Int                                       # fixed strength, or index into FPlant.ε (0 = fixed)
end
struct LineSpec
    ix::Int; other::Int; α::Float64; logistic::Bool; r::Float64; K::Float64; θ::Float64; d::Float64
    mods::Vector{ModSpec}
end
function line_specs(base::Plant, iε)
    S = base.sys; p = base.p; out = LineSpec[]
    for (j, l) in enumerate(S.envs[1].lines)
        L = S.lines[l]; xi = S.xi[1][j]; off = S.hb[1][j]
        (hasproperty(L, :chi) && sv(L.chi, p, Float64) != 0) && error("χ interaction not supported in the fast plant")
        L.growth in (:logistic, :theta) || error("growth law $(L.growth) not supported in the fast plant")
        mods = ModSpec[]
        for M in modules(L)
            ncell(M.action) == 0 && naux(M.action) == 0 || error("action $(M.action) not supported in the fast plant")
            M.kernel in (:none, :chain3) || error("kernel $(M.kernel) not supported in the fast plant")
            M.action in (:hill, :capacity, :loss) || error("action $(M.action) not supported in the fast plant")
            nk = nker(M.kernel); gs(k) = hasproperty(M.s, k) ? sv(getfield(M.s, k), p, Float64) : 0.0
            εs = M.s.ε; εj = εs.i > 0 ? findfirst(==(εs.i), iε) : 0
            εs.i > 0 && εj === nothing && error("unexpected free slot")
            push!(mods, ModSpec(M.input == :D, nk, nk > 0 ? gs(:k) : 0.0, off, M.action, gs(:β), gs(:h), gs(:c), gs(:d0), ec50_of(M.s, p, Float64),
                                εs.i > 0 ? NaN : εs.v, εj))
            off += mwidth(M)
        end
        m = 3 - j
        push!(out, LineSpec(xi + 1, S.xi[1][m] + 1, sv(S.alpha[l][S.envs[1].lines[m]], p, Float64), L.growth == :logistic,
                            sv(L.s.r, p, Float64), sv(L.s.K, p, Float64), sv(L.s.θ, p, Float64), sv(L.s.d, p, Float64), mods))
    end
    out
end

struct FPlant{TT<:NamedTuple,HT<:NamedTuple}
    name::String
    base::Plant                 # sys with the low-resource modules; p = [scenario params; ε_N, ε_C]
    iε::Vector{Int}             # positions of ε_N, ε_C in base.p
    ε::Vector{Float64}          # Stage 4 (or bootstrap) ε_N, ε_C
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
    c4 = zgroup("s4").cands[1]
    ε = b === nothing ? copy(c4.params) : ZOO_STATE.results["s4"].res.pars[1][b, :]
    any(isnan, ε) && return nothing
    np = length(P0.p)
    lines = map(1:2) do j
        L4 = c4.sys.lines[j]
        merge(P0.sys.lines[j], (extra=[merge(M, (s=shift_slots(M.s, np),)) for M in L4.extra],))
    end
    p = vcat(P0.p, ε)
    base = make_plant(P0.name, lines, P0.sys.alpha, p; uR=0.0, tox=tox)
    n = base.sys.n
    name = (scen == :A ? "A: Stage 3 coculture model as selected" : "B: monoculture drug response for both lines") *
           " + fitted low-resource modules"
    FPlant(name, base, [np + 1, np + 2], ε, eps_mode, Float64(FCFG["eps_u_ref"]),
           n + 4, base.iN, base.iC, n + 1, n + 2, n + 3, n + 4, base.K, tox, host, line_specs(base, [np + 1, np + 2]), Float64(koff))
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
        crowd = max(x[xi+1], zero(T))
        for (m, lm) in enumerate(env.lines)
            m == j && continue
            crowd += sv(S.alpha[l][lm], p, T) * max(x[S.xi[1][m]+1], zero(T))
        end
        line_rhs_multi!(du, x, xi, S.hb[1][j], L, p, envc, crowd, T)
    end
    Bd = x[FP.iB]; Tx = x[FP.iT]; Pr = x[FP.iP]; Wt = x[FP.iW]; tp = FP.tox; h = FP.host
    du[FP.iB] = tp.s_B * u - tp.d_B * Bd
    du[FP.iT] = tp.a_H * (1 - h.p_H * Pr) * (1 - exp(-max(Bd, zero(T)))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
    du[FP.iP] = h.k_P * (v - Pr)
    du[FP.iW] = -h.k_W * v + h.r_W * (1 - Wt) * (1 - v)
    nothing
end
function fast_rhs!(du, x, w, FP::FPlant)
    u, v = w[1], w[2]
    T = promote_type(eltype(x), eltype(w)); z = zero(T)
    f = εfactor(FP, T(u)); uR = T(U_R) * v
    @inbounds for L in FP.spec
        X = max(x[L.ix], z); crowd = X + L.α * max(x[L.other], z)
        Keff = T(L.K); loss = T(L.d)
        for M in L.mods
            ε = M.εj == 0 ? T(M.ε) : 1 + (FP.ε[M.εj] - 1) * f
            wi = (M.input_D ? T(u) : uR) * ε
            E = wi
            if M.kernel_chain == 3
                o = M.o
                # low-resource stress chains may fade faster after refeeding (sensitivity only; default = fitted rate)
                k = (M.input_D || isnan(FP.koff)) ? T(M.k) : M.k * v + FP.koff * (1 - v)
                du[o+1] = k * (wi - x[o+1]); du[o+2] = k * (x[o+1] - x[o+2]); du[o+3] = k * (x[o+2] - x[o+3])
                E = max(x[o+3], z)
            end
            if M.action == :hill
                loss += E > 0 ? M.β * E^M.h / (T(M.ec50)^M.h + E^M.h) : z
            elseif M.action == :capacity
                Keff /= 1 + M.c * E
            else                                     # :loss
                loss += M.d0 * wi + M.β * E
            end
        end
        q = max(crowd / Keff, z)
        du[L.ix] = L.r * X * (L.logistic ? 1 - q : 1 - q^L.θ) - loss * X
    end
    Bd = x[FP.iB]; Tx = x[FP.iT]; Pr = x[FP.iP]; Wt = x[FP.iW]; tp = FP.tox; h = FP.host
    du[FP.iB] = tp.s_B * u - tp.d_B * Bd
    du[FP.iT] = tp.a_H * (1 - h.p_H * Pr) * (1 - exp(-max(Bd, z))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
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
function linearize2(x, w, FP::FPlant, Δ)
    nx = length(x)
    J = ForwardDiff.jacobian(z -> F2(z[1:nx], z[nx+1:end], FP, Δ), vcat(x, w))
    F2(x, w, FP, Δ), J[:, 1:nx], J[:, nx+1:end]
end
function x0_fast(FP::FPlant, N0, C0)
    x = zeros(FP.nx); x[FP.iN] = N0; x[FP.iC] = C0; x[FP.iW] = 1.0
    x
end

# ── cost ─────────────────────────────────────────────────────────────────────
struct Weights2
    qN::Float64; qC::Float64; qT::Float64; qW::Float64; ru::Float64; rv::Float64; tf::Float64
    μC::Float64; μT::Float64; μW::Float64; δ::Float64
    Bc::Float64; Btox::Float64; Wmin::Float64; umax::Float64; K::Float64
end
function weights2(FP::FPlant; umax=Float64(CFG["dose"]["u_max_data"]), Wmin=Float64(HCFG["W_min"]), w=nothing, rv=Float64(FCFG["rv"]))
    c = CFG["weights"]
    qN, qC, qT, ru = w === nothing ? (c["qN"], c["qC"], c["qT"], c["ru"]) : Tuple(w)
    Weights2(qN, qC, qT, HCFG["qW"], ru, rv, c["terminal_factor"], c["mu_cancer"], c["mu_tox"], HCFG["mu_host"], c["barrier_delta"],
             CFG["cancer"]["B_cancer_fraction_of_K"] * FP.K, CFG["toxicity"]["B_tox"], Wmin, umax, FP.K)
end
wloss(x, FP, W) = (1 - x[FP.iW]) / (1 - W.Wmin)
function stage_cost2(x, w, FP::FPlant, W::Weights2)
    N, C, T = x[FP.iN], x[FP.iC], x[FP.iT]
    W.qN * (N / W.K)^2 + W.qC * (C / W.K)^2 + W.qT * (T / W.Btox)^2 + W.qW * wloss(x, FP, W)^2 +
        W.ru * (w[1] / W.umax)^2 + W.rv * w[2]^2 +
        relaxed_barrier(1 - (N + C) / W.Bc, W.μC, W.δ) + relaxed_barrier(1 - T / W.Btox, W.μT, W.δ) +
        relaxed_barrier(1 - wloss(x, FP, W), W.μW, W.δ)
end
terminal_cost2(x, FP::FPlant, W::Weights2) = W.tf * (W.qN * (x[FP.iN] / W.K)^2 + W.qC * (x[FP.iC] / W.K)^2 + W.qT * (x[FP.iT] / W.Btox)^2 + W.qW * wloss(x, FP, W)^2)

function total_cost2(x0, U, FP, Δ, W)
    x = copy(x0); J = 0.0
    for k in axes(U, 1)
        x = F2(x, U[k, :], FP, Δ); J += stage_cost2(x, U[k, :], FP, W)
    end
    J + terminal_cost2(x, FP, W)
end

# adjoint (vector–Jacobian) gradient and Gauss–Newton Hessian of J over the stacked controls
function grad_hess2(x0, U, FP::FPlant, Δ, W::Weights2)
    H = size(U, 1); m = size(U, 2)
    xs = [copy(x0)]; As = Matrix{Float64}[]; Bs = Matrix{Float64}[]
    for k in 1:H
        xn, A, B = linearize2(xs[k], U[k, :], FP, Δ)
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
function discrete_solve2(x0, FP::FPlant, Δ, H, W::Weights2; ulevels=Float64.(CFG["dose"]["levels_data"]),
                         vlevels=Float64.(FCFG["levels"]), maxrun=Int(FCFG["max_consecutive_intervals"]),
                         mingap=Int(FCFG["min_gap_intervals"]), run0=0, gap0=10^6, Uinit=nothing)
    acts = [(u, v) for v in vlevels for u in ulevels]; na = length(acts); nx = length(x0)
    buf = RK4Buf(nx)
    X = [[zeros(nx) for _ in 1:na] for _ in 1:H]; cst = [zeros(na) for _ in 1:H]; ord = [collect(1:na) for _ in 1:H]
    best = Ref(Inf); bestU = zeros(H, 2); U = zeros(H, 2)
    if Uinit !== nothing                      # upper bound from the warm-start plan, if it obeys the rules
        x = copy(x0); c = 0.0; run, gap = run0, gap0; ok = true; y = zeros(nx)
        for k in 1:H
            u, v = Uinit[k, 1], Uinit[k, 2]
            admissible(v, run, gap, maxrun, mingap) || (ok = false; break)
            F2!(y, x, (u, v), FP, Δ, buf); x, y = y, x; c += stage_cost2(x, (u, v), FP, W)
            run, gap = v > 0 ? (run + 1, 0) : (0, gap + 1)
        end
        ok && (best[] = c + terminal_cost2(x, FP, W); bestU .= Uinit)
    end
    function dfs(x, k, cost, run, gap)
        if k > H
            c = cost + terminal_cost2(x, FP, W)
            c < best[] && (best[] = c; bestU .= U)
            return
        end
        for a in 1:na
            u, v = acts[a]
            if admissible(v, run, gap, maxrun, mingap)
                F2!(X[k][a], x, (u, v), FP, Δ, buf); cst[k][a] = cost + stage_cost2(X[k][a], (u, v), FP, W)
            else
                cst[k][a] = Inf
            end
        end
        sortperm!(ord[k], cst[k])
        for a in ord[k]
            cst[k][a] >= best[] && break
            u, v = acts[a]; U[k, 1] = u; U[k, 2] = v
            dfs(X[k][a], k + 1, cst[k][a], v > 0 ? run + 1 : 0, v > 0 ? 0 : gap + 1)
        end
    end
    dfs(copy(x0), 1, 0.0, run0, gap0)
    (U=copy(bestU), J=best[])
end

# receding horizon on plant Pp, planning with Pm; the fasting rules carry over between decisions
function closed_loop2(Pm::FPlant, Pp::FPlant, x0, Δ, n_dec, H, W::Weights2; mode=:discrete, shrink=true, kw...)
    x = copy(x0); applied = zeros(n_dec, 2); run = 0; gap = 10^6
    Uwarm = repeat([mode == :discrete ? 0.0 : 0.5  0.0], H)
    for k in 1:n_dec
        Hk = shrink ? min(H, n_dec - k + 1) : H
        sol = mode == :discrete ? discrete_solve2(x, Pm, Δ, Hk, W; run0=run, gap0=gap, Uinit=Uwarm[1:Hk, :], kw...) :
                                  sqp_solve2(x, Uwarm[1:Hk, :], Pm, Δ, W)
        w = sol.U[1, :]; applied[k, :] = w
        run = w[2] > 0 ? run + 1 : 0; gap = w[2] > 0 ? 0 : gap + 1
        x = F2(x, w, Pp, Δ)
        Uwarm = vcat(sol.U[2:end, :], repeat([sol.U[end, 1] 0.0], H - Hk + 1))[1:H, :]     # shift; pad fed
    end
    (U=applied,)
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
]

function outcomes2(FP::FPlant, x0, U, Δ, W::Weights2)
    ts, xs = trajectory2(x0, U, FP, Δ)
    N = [x[FP.iN] for x in xs]; C = [x[FP.iC] for x in xs]; T = [x[FP.iT] for x in xs]
    Pr = [x[FP.iP] for x in xs]; Wt = [x[FP.iW] for x in xs]; tot = N .+ C
    over = findfirst(>=(W.Bc), tot)
    (N_final=N[end], C_final=C[end], total_final=tot[end], resistant_fraction=C[end] / max(tot[end], 1e-9),
     max_total=maximum(tot), max_tox=maximum(T), min_weight=minimum(Wt), cum_dose=sum(U[:, 1]) * Δ, fast_days=sum(U[:, 2]) * Δ,
     cancer_violation=over !== nothing, tox_violation=maximum(T) >= W.Btox, weight_violation=minimum(Wt) <= W.Wmin,
     ts=ts, N=N, C=C, T=T, P=Pr, W=Wt)
end
