# Sequential (successive-linearisation) model-predictive control of cisplatin dosing.
# Include after control_model.jl.
#
# Finite-horizon problem at every decision time (decision interval Δ, H intervals):
#   minimise  J(U) = Σ_{k=1..H} ℓ(x_{k+1}, u_k) + φ(x_{H+1}),   x_{k+1} = F(x_k, u_k)
#   subject to 0 ≤ u_k ≤ u_max  [and Δ·Σu_k ≤ U_cum]
#   ℓ = qB(N+C)/K + qR·C/K + qT(T/B_tox)² + ru(u/u_max)² + b_μC(1 − (N+C)/B_cancer) + b_μT(1 − T/B_tox)
#   φ = terminal_factor · [qB(N+C)/K + qR·C/K + qT(T/B_tox)²]
# Goal (user decision 2026-10-02): low total burden and resistance held back. The burden is linear in
# the counts (every cell costs the same at any tumour size); qR is the extra cost of a resistant cell.
# b_μ is a relaxed logarithmic barrier (Feller & Ebenbauer 2017, IEEE Trans Autom Control 62:1223):
#   −μ log z for z ≥ δ, and its quadratic continuation below δ, so an infeasible start still has a
#   finite, steep cost and the optimiser is pushed back inside instead of failing.
#
# Continuous mode ("demonstration"): SQP iterations. Each iteration
#   1. rolls the nonlinear model forward,
#   2. linearises F along the rollout (A_k, B_k by forward-mode AD, control_model.jl: linearize),
#   3. computes ∇J by the adjoint recursion μ_k = ∇ℓ_k + A_kᵀ μ_{k+1} (vector–Jacobian products),
#   4. builds a Gauss–Newton Hessian from the forward sensitivities S_{k+1} = A_k S_k + B_k e_kᵀ and
#      the PSD part of ∇²ℓ (ForwardDiff Hessians),
#   5. solves the box-constrained QP for the step (projected accelerated gradient; the dose budget,
#      if any, by bisection on its multiplier),
#   6. line-searches on the true nonlinear cost.
# This is SQP / sequential linearisation MPC, not a fixed LQR: the local QP has the LQR structure
# (δx_{k+1} = A_k δx_k + B_k δu_k with quadratic cost) but is re-linearised every iteration and
# every decision time, has bounds on u, and carries non-quadratic barrier costs.
#
# Discrete mode ("data-supported"): u restricted to the measured doses; the finite-horizon problem is
# solved exactly by depth-first search over all dose sequences, sharing rollout prefixes and pruning
# branches whose partial cost already exceeds the best (every stage cost is ≥ 0).

struct Weights
    qB::Float64; qR::Float64; qT::Float64; ru::Float64; tf::Float64
    μC::Float64; μT::Float64; δ::Float64
    Bc::Float64; Btox::Float64; umax::Float64; K::Float64
end
function default_weights(P::Plant; umax=Float64(CFG["dose"]["u_max_data"]), w=nothing, Btox=Float64(CFG["toxicity"]["B_tox"]))
    c = CFG["weights"]
    qB, qR, qT, ru = w === nothing ? (c["qB"], c["qR"], c["qT"], c["ru"]) : Tuple(w)
    Weights(qB, qR, qT, ru, c["terminal_factor"], c["mu_cancer"], c["mu_tox"], c["barrier_delta"],
            CFG["cancer"]["B_cancer_fraction_of_K"] * P.K, Btox, umax, P.K)
end

# tumour burden plus the resistance premium, linear in the counts (≥ 0, so the exact search can prune)
burden_cost(N, C, W) = (W.qB * (max(N, zero(N)) + max(C, zero(C))) + W.qR * max(C, zero(C))) / W.K
relaxed_barrier(z, μ, δ) = z >= δ ? -μ * log(z) : μ * (0.5 * ((z - 2δ) / δ)^2 - 0.5 - log(δ))
function stage_cost(x, u, P::Plant, W::Weights)
    N, C, T = x[P.iN], x[P.iC], x[P.iT]
    burden_cost(N, C, W) + W.qT * (T / W.Btox)^2 + W.ru * (u / W.umax)^2 +
        relaxed_barrier(1 - (N + C) / W.Bc, W.μC, W.δ) + relaxed_barrier(1 - T / W.Btox, W.μT, W.δ)
end
terminal_cost(x, P::Plant, W::Weights) = W.tf * (burden_cost(x[P.iN], x[P.iC], W) + W.qT * (x[P.iT] / W.Btox)^2)

function rollout(x0, U, P::Plant, Δ)
    xs = [copy(x0)]
    for u in U
        push!(xs, F(xs[end], u, P, Δ))
    end
    xs
end
function total_cost(x0, U, P::Plant, Δ, W::Weights)
    xs = rollout(x0, U, P, Δ)
    sum(stage_cost(xs[k+1], U[k], P, W) for k in eachindex(U)) + terminal_cost(xs[end], P, W)
end

# gradient of J by the adjoint (reverse) recursion over the explicit Jacobians: vector–Jacobian products
function adjoint_gradient(xs, As, Bs, U, P, W)
    H = length(U)
    gx(k) = ForwardDiff.gradient(x -> stage_cost(x, U[k], P, W), xs[k+1])
    gu(k) = ForwardDiff.derivative(u -> stage_cost(xs[k+1], u, P, W), U[k])
    μ = gx(H) .+ ForwardDiff.gradient(x -> terminal_cost(x, P, W), xs[H+1])     # ∂J/∂x_{H+1}
    g = zeros(H)
    for k in H:-1:1
        g[k] = gu(k) + LA.dot(Bs[k], μ)                                         # ∂J/∂u_k = ℓ_u + B_kᵀ μ_{k+1}
        k > 1 && (μ = gx(k - 1) .+ As[k]' * μ)                                  # μ_k = ∇ℓ(x_k) + A_kᵀ μ_{k+1}
    end
    g
end

psd(M) = (E = LA.eigen(LA.Symmetric((M + M') / 2)); E.vectors * LA.Diagonal(max.(E.values, 0.0)) * E.vectors')

# Gauss–Newton Hessian of J from forward sensitivities
function gn_hessian(xs, As, Bs, U, P, W)
    H = length(U); nx = length(xs[1])
    S = zeros(nx, H); Hs = zeros(H, H)
    for k in 1:H
        S = As[k] * S; S[:, k] .+= Bs[k]                                        # S_{k+1} = A_k S_k + B_k e_kᵀ
        Q = ForwardDiff.hessian(x -> stage_cost(x, U[k], P, W), xs[k+1])
        k == H && (Q = Q + ForwardDiff.hessian(x -> terminal_cost(x, P, W), xs[k+1]))
        Hs .+= S' * psd(Q) * S
        Hs[k, k] += max(ForwardDiff.derivative(u -> ForwardDiff.derivative(v -> stage_cost(xs[k+1], v, P, W), u), U[k]), 0.0)
    end
    (Hs + Hs') / 2 + 1e-9 * LA.I
end

# box QP  min ½vᵀHv + cᵀv, lo ≤ v ≤ hi  (projected accelerated gradient)
function box_qp(Hm, c, lo, hi, v0; iters=3000)
    L = max(LA.opnorm(Hm), 1e-12); v = clamp.(v0, lo, hi); y = copy(v); t = 1.0
    for _ in 1:iters
        vn = clamp.(y .- (Hm * y .+ c) ./ L, lo, hi)
        tn = (1 + sqrt(1 + 4t^2)) / 2
        y = vn .+ ((t - 1) / tn) .* (vn .- v)
        maximum(abs.(vn .- v)) < 1e-12 && (v = vn; break)
        v = vn; t = tn
    end
    v
end
# with an optional budget Δ·Σv ≤ cum (bisection on its multiplier)
function qp_step(Hm, g, U, lo, hi, Δ, cum)
    c = g .- Hm * U
    v = box_qp(Hm, c, lo, hi, U)
    (cum <= 0 || Δ * sum(v) <= cum + 1e-12) && return v
    a, b = 0.0, 1.0
    while Δ * sum(box_qp(Hm, c .+ b * Δ, lo, hi, U)) > cum; b *= 2; end
    for _ in 1:60
        m = (a + b) / 2
        Δ * sum(box_qp(Hm, c .+ m * Δ, lo, hi, U)) > cum ? (a = m) : (b = m)
    end
    box_qp(Hm, c .+ b * Δ, lo, hi, U)
end

# SQP solve of the finite-horizon problem from x0 (continuous doses)
function sqp_solve(x0, U0, P::Plant, Δ, W::Weights; cum=0.0, max_iter=Int(CFG["sqp"]["max_iter"]), tol=Float64(CFG["sqp"]["tol"]))
    H = length(U0); lo = zeros(H); hi = fill(W.umax, H)
    U = clamp.(copy(U0), lo, hi)
    if cum > 0 && Δ * sum(U) > cum
        U .*= cum / (Δ * sum(U))
    end
    J = total_cost(x0, U, P, Δ, W); hist = [J]; it = 0; Jold = J
    for outer it in 1:max_iter
        xs = [copy(x0)]; As = Matrix{Float64}[]; Bs = Vector{Float64}[]
        for k in 1:H
            xn, A, B = linearize(xs[k], U[k], P, Δ)
            push!(xs, xn); push!(As, A); push!(Bs, B)
        end
        g = adjoint_gradient(xs, As, Bs, U, P, W)
        Hm = gn_hessian(xs, As, Bs, U, P, W)
        d = qp_step(Hm, g, U, lo, hi, Δ, cum) .- U
        maximum(abs.(d)) < 1e-8 && break
        α = 1.0; slope = LA.dot(g, d); accepted = false
        for _ in 1:25
            Un = U .+ α .* d; Jn = total_cost(x0, Un, P, Δ, W)
            if Jn <= J + 1e-4 * α * min(slope, 0.0)
                U, Jold, J = Un, J, Jn; accepted = true; break
            end
            α /= 2
        end
        push!(hist, J)
        (!accepted || abs(Jold - J) < tol * (1 + abs(J))) && break
    end
    (U=U, J=J, iters=it, history=hist)
end

# exact search over discrete dose levels (prefix sharing + pruning)
function discrete_solve(x0, P::Plant, Δ, H, levels, W::Weights; min_gap=0.0, cum=0.0)
    best = Ref(Inf); bestU = Ref(zeros(H))
    gap_steps = ceil(Int, min_gap / Δ - 1e-9)
    function dfs(x, k, U, cost, last_dose_k, used)
        if k > H
            c = cost + terminal_cost(x, P, W)
            c < best[] && (best[] = c; bestU[] = copy(U))
            return
        end
        for u in levels
            u > 0 && gap_steps > 0 && k - last_dose_k <= gap_steps && continue
            (cum > 0 && used + u * Δ > cum + 1e-12) && continue
            xn = F(x, u, P, Δ)
            c = cost + stage_cost(xn, u, P, W)
            c >= best[] && continue
            U[k] = u
            dfs(xn, k + 1, U, c, u > 0 ? k : last_dose_k, used + u * Δ)
        end
    end
    dfs(copy(x0), 1, zeros(H), 0.0, -10^6, 0.0)
    (U=bestU[], J=best[])
end

# receding-horizon closed loop on the plant Pp, planning with Pm (the same model unless testing mismatch)
function closed_loop(Pm::Plant, Pp::Plant, x0, Δ, n_dec, H, W::Weights; mode=:continuous, levels=Float64[],
                     shrink=true, cum=0.0, min_gap=0.0, u_init=0.5)
    x = copy(x0); applied = Float64[]; Uwarm = fill(min(u_init, W.umax), H); used = 0.0; iters = Int[]
    for k in 1:n_dec
        Hk = shrink ? min(H, n_dec - k + 1) : H
        rem = cum > 0 ? max(cum - used, 0.0) : 0.0
        if mode == :continuous
            sol = sqp_solve(x, Uwarm[1:Hk], Pm, Δ, W; cum=cum > 0 ? max(rem, 1e-9) : 0.0)
            push!(iters, sol.iters)
        else
            sol = discrete_solve(x, Pm, Δ, Hk, levels, W; min_gap=min_gap, cum=cum > 0 ? rem : 0.0)
        end
        u = sol.U[1]; push!(applied, u); used += u * Δ
        x = F(x, u, Pp, Δ)
        Uwarm = vcat(sol.U[2:end], fill(sol.U[end], H - Hk + 1))[1:H]
    end
    (U=applied, iters=iters)
end

# ── baseline policies (state feedback at decision times) ────────────────────
function run_policy(policy, P::Plant, x0, Δ, n_dec)
    x = copy(x0); U = Float64[]; on = true
    for k in 1:n_dec
        u = policy(x, k)
        push!(U, u); x = F(x, u, P, Δ)
    end
    U
end
# Threshold adaptive therapy in the spirit of Zhang et al. 2017 (Nat Commun 8:1816): treat until the
# burden falls to a lower threshold, pause until it regrows to an upper one. Thresholds here are set
# relative to B_cancer because the cultures start from seeding, not from an established tumour.
function threshold_policy(P::Plant, W::Weights; dose=1.0, hi=0.5, lo=0.25)
    on = Ref(false)
    (x, k) -> begin
        tot = x[P.iN] + x[P.iC]
        tot >= hi * W.Bc && (on[] = true)
        tot <= lo * W.Bc && (on[] = false)
        on[] ? dose : 0.0
    end
end
baseline_policies(P::Plant, W::Weights) = [
    ("No treatment", (x, k) -> 0.0),
    ("Constant 1 µM", (x, k) -> 1.0),
    ("Maximum measured dose (1.47 µM)", (x, k) -> 1.47),
    ("Periodic 1 µM, 2 days on / 2 off", (x, k) -> isodd(k) ? 1.0 : 0.0),
    ("Threshold adaptive (Zhang-style)", threshold_policy(P, W)),
]

# biologically interpretable outcomes of a schedule
function outcomes(P::Plant, x0, U, Δ, W::Weights)
    ts, xs = trajectory(x0, U, P, Δ)
    N = [x[P.iN] for x in xs]; C = [x[P.iC] for x in xs]; T = [x[P.iT] for x in xs]; tot = N .+ C
    over = findfirst(>=(W.Bc), tot)
    (N_final=N[end], C_final=C[end], total_final=tot[end], resistant_fraction=C[end] / max(tot[end], 1e-9),
     max_total=maximum(tot), max_tox=maximum(T), cum_dose=sum(U) * Δ, cancer_violation=over !== nothing,
     tox_violation=maximum(T) >= W.Btox, day_over_bound=over === nothing ? NaN : ts[over], ts=ts, N=N, C=C, T=T)
end
