# Sequential (successive-linearisation) model-predictive control of cisplatin dosing.
# Include after control_model.jl.
#
# Finite-horizon problem at every decision time (decision interval Δ, H intervals):
#   minimise  J(U) = Σ_{k=1..H} ℓ(x_{k+1}, u_k) + φ(x_{H+1}),   x_{k+1} = F(x_k, u_k)
#   subject to 0 ≤ u_k ≤ u_max  [and Δ·Σu_k ≤ U_cum]
#   ℓ = qP·σ((N+C − B_prog)/s) + qB(N+C)/K + qR·C/K + qT(T/B_tox)² + ru(u/u_max)²
#       + b_μC(1 − (N+C)/B_cancer) + b_μT(1 − T/B_tox)
#   φ = terminal_factor · [qP·σ(…) + qB(N+C)/K + qR·C/K + qT(T/B_tox)²] + qRT·C/K
# Goal (user decision 2026-10-07): TIME TO PROGRESSION. σ is a smooth step (logistic, width s) at the
# progression threshold B_prog (80% of K), so Σℓ counts the decision intervals spent progressed: the
# controller keeps the burden under B_prog for as long as it can, e.g. treating hard when the culture
# nears 80% confluence and then holding off. qRT·C/K at the end of the look-ahead prices the resistant
# cells left behind, the cells that will cause the next, untreatable progression. The burden and
# resistance terms of the 2026-10-02 objective (qB, qR) remain as small tie-breakers (swept).
# b_μ is a relaxed logarithmic barrier (Feller & Ebenbauer 2017, IEEE Trans Autom Control 62:1223):
#   −μ log z for z ≥ δ, and its quadratic continuation below δ, so a state beyond the limit has a
#   finite, steep cost. It is a SOFT constraint (a penalty), not a hard limit: at B_cancer it is the
#   "soft burden penalty", which every growing normal-medium culture exceeds (counted in the report).
#   B_cancer_fraction_of_K = 0 (sweep only) removes it: Bc = Inf, so b_μC = −μ log 1 = 0.
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
# branches whose partial cost already exceeds the best (every stage cost is ≥ 0). The search starts
# from the cost of the warm-start plan (the previous plan shifted by one interval) and expands the
# children of each node cheapest-first, so it prunes early; it stays exact.
#
# Receding horizon: the look-ahead is FIXED, also for the last decisions of a run (the plant is simply
# simulated past the run's end), so a run's schedule is the first decisions of a continuing policy and
# does not depend on the run length (no end-of-run effect).
# Move blocking (Cagienard et al. 2007, Automatica 43:1121): the look-ahead is a list of blocks of
# decision intervals with one dose per block (control_config.toml lookahead_blocks, e.g. [1,1,2,2,4,4] =
# 28 days in 6 decisions), so a long look-ahead keeps the exact search affordable. The first block is
# one interval: the applied dose is always decided at full resolution. Continuous (SQP) mode plans every
# interval of the same look-ahead.

import Random

struct Weights
    qB::Float64; qR::Float64; qT::Float64; ru::Float64; tf::Float64
    μC::Float64; μT::Float64; δ::Float64
    Bc::Float64; Btox::Float64; umax::Float64; K::Float64
    qP::Float64; Bp::Float64; sP::Float64; qRT::Float64      # progression weight, threshold, step width; terminal resistant weight
end
# a copy of a weight struct with some fields replaced, e.g. with_fields(W; Bc=…)
with_fields(W::T; kw...) where {T} = T((haskey(kw, f) ? kw[f] : getfield(W, f) for f in fieldnames(T))...)
# the objective weights (qB, qR, qT, ru, qP, qRT): w overrides them (a 4-vector keeps the default qP, qRT)
function objective_weights(w)
    c = CFG["weights"]; p = CFG["progression"]
    w === nothing && return (c["qB"], c["qR"], c["qT"], c["ru"], p["qP"], p["qRT"])
    length(w) == 6 ? Tuple(Float64.(w)) : (Float64.(w)..., p["qP"], p["qRT"])
end
# progression threshold B_prog and step width s (cells) for capacity K; prog_frac = 0 removes the term
function progression_scale(K; prog_frac=Float64(CFG["progression"]["threshold_fraction_of_K"]))
    (prog_frac > 0 ? prog_frac * K : Inf, Float64(CFG["progression"]["step_width_fraction_of_K"]) * K)
end
# bc_frac: B_cancer as a fraction of K (0 = no burden penalty); tf: terminal factor (both swept in the report)
function default_weights(P::Plant; umax=Float64(CFG["dose"]["u_max_data"]), w=nothing, Btox=Float64(CFG["toxicity"]["B_tox"]),
                         bc_frac=Float64(CFG["cancer"]["B_cancer_fraction_of_K"]), tf=Float64(CFG["weights"]["terminal_factor"]),
                         prog_frac=Float64(CFG["progression"]["threshold_fraction_of_K"]))
    c = CFG["weights"]
    qB, qR, qT, ru, qP, qRT = objective_weights(w)
    Bp, sP = progression_scale(P.K; prog_frac=prog_frac)
    Weights(qB, qR, qT, ru, tf, c["mu_cancer"], c["mu_tox"], c["barrier_delta"],
            bc_frac > 0 ? bc_frac * P.K : Inf, Btox, umax, P.K, qP, Bp, sP, qRT)
end

# tumour burden plus the resistance premium, linear in the counts (≥ 0, so the exact search can prune)
burden_cost(N, C, W) = (W.qB * (max(N, zero(N)) + max(C, zero(C))) + W.qR * max(C, zero(C))) / W.K
relaxed_barrier(z, μ, δ) = z >= δ ? -μ * log(z) : μ * (0.5 * ((z - 2δ) / δ)^2 - 0.5 - log(δ))
# smooth step: ≈ 0 below the progression threshold, ≈ 1 above it (logistic of width s, overflow-safe)
function progressed(N, C, W)
    z = (max(N, zero(N)) + max(C, zero(C)) - W.Bp) / W.sP
    z >= 0 ? 1 / (1 + exp(-z)) : exp(z) / (1 + exp(z))
end
progression_cost(N, C, W) = W.qP * progressed(N, C, W)
function stage_cost(x, u, P::Plant, W::Weights)
    N, C, T = x[P.iN], x[P.iC], x[P.iT]
    progression_cost(N, C, W) + burden_cost(N, C, W) + W.qT * (T / W.Btox)^2 + W.ru * (u / W.umax)^2 +
        relaxed_barrier(1 - (N + C) / W.Bc, W.μC, W.δ) + relaxed_barrier(1 - T / W.Btox, W.μT, W.δ)
end
# terminal: the state-cost terms ×tf, plus the resistant cells left at the end of the look-ahead
terminal_cost(x, P::Plant, W::Weights) = W.tf * (progression_cost(x[P.iN], x[P.iC], W) + burden_cost(x[P.iN], x[P.iC], W) + W.qT * (x[P.iT] / W.Btox)^2) +
    W.qRT * max(x[P.iC], zero(x[P.iC])) / W.K

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

# The stage cost ℓ split into its terms, summed over the intervals of an applied schedule (the terminal
# cost is a planning device and is left out). Same terms and weights as stage_cost: total = Σ_k ℓ.
function cost_terms(P::Plant, x0, U, Δ, W::Weights)
    xs = rollout(x0, U, P, Δ); z = zeros(7)
    for k in eachindex(U)
        x = xs[k+1]; N, C, T = x[P.iN], x[P.iC], x[P.iT]
        z .+= [W.qB * (max(N, 0.0) + max(C, 0.0)) / W.K, W.qR * max(C, 0.0) / W.K, W.qT * (T / W.Btox)^2, W.ru * (U[k] / W.umax)^2,
               relaxed_barrier(1 - (N + C) / W.Bc, W.μC, W.δ), relaxed_barrier(1 - T / W.Btox, W.μT, W.δ), progression_cost(N, C, W)]
    end
    (burden=z[1], resistance=z[2], toxicity=z[3], dose=z[4], penalty_burden=z[5], barrier_tox=z[6], progression=z[7], total=sum(z))
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
# ∇J at a schedule U from x0 (rollout, linearisation, adjoint); used by the optimality test
function cost_gradient(x0, U, P::Plant, Δ, W::Weights)
    xs = [copy(x0)]; As = Matrix{Float64}[]; Bs = Vector{Float64}[]
    for k in eachindex(U)
        xn, A, B = linearize(xs[k], U[k], P, Δ); push!(xs, xn); push!(As, A); push!(Bs, B)
    end
    adjoint_gradient(xs, As, Bs, U, P, W)
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

# look-ahead blocks: an Int H means H single intervals; a vector lists the intervals of each decision block
as_blocks(H::Integer) = fill(1, H)
as_blocks(H::AbstractVector) = collect(Int, H)
# the blocks cut to at most m intervals in total (the old shrinking horizon; comparison only)
function cut_blocks(b, m)
    out = Int[]; left = m
    for x in b
        left <= 0 && break
        push!(out, min(x, left)); left -= x
    end
    out
end
expand_blocks(Ub, b) = reduce(vcat, [fill(Ub[j], b[j]) for j in eachindex(b)])   # one dose per interval

# exact search over discrete dose levels (prefix sharing, warm-start bound, cheapest-first, pruning).
# H: an Int (H single intervals) or the look-ahead blocks; a block applies one dose to its intervals, and
# the cost is summed over every interval, so blocking changes the plan's resolution, not the objective.
function discrete_solve(x0, P::Plant, Δ, H, levels, W::Weights; min_gap=0.0, cum=0.0, Uinit=nothing)
    blk = as_blocks(H); H = length(blk)
    (all(==(1), blk) || (min_gap <= 0 && cum <= 0)) || error("dose gaps and budgets need single-interval blocks")
    step(x, u, k) = (c = 0.0; for _ in 1:blk[k]; x = F(x, u, P, Δ); c += stage_cost(x, u, P, W); end; (x, c))
    na = length(levels); best = Ref(Inf); bestU = zeros(H); U = zeros(H)
    gap_steps = ceil(Int, min_gap / Δ - 1e-9)
    allowed(u, k, last, used) = !(u > 0 && gap_steps > 0 && k - last <= gap_steps) && !(cum > 0 && used + u * Δ > cum + 1e-12)
    if Uinit !== nothing && length(Uinit) >= H           # upper bound from the warm-start plan, if it is admissible
        x = copy(x0); c = 0.0; last = -10^6; used = 0.0; ok = true
        for k in 1:H
            u = Uinit[k]
            (any(==(u), levels) && allowed(u, k, last, used)) || (ok = false; break)
            x, ck = step(x, u, k); c += ck
            u > 0 && (last = k); used += u * Δ
        end
        ok && (best[] = c + terminal_cost(x, P, W); bestU .= Uinit[1:H])
    end
    X = [Vector{Vector{Float64}}(undef, na) for _ in 1:H]; cst = [zeros(na) for _ in 1:H]; ord = [collect(1:na) for _ in 1:H]
    function dfs(x, k, cost, last, used)
        if k > H
            c = cost + terminal_cost(x, P, W)
            c < best[] && (best[] = c; bestU .= U)
            return
        end
        for a in 1:na
            u = levels[a]
            if allowed(u, k, last, used)
                X[k][a], ck = step(x, u, k); cst[k][a] = cost + ck
            else
                cst[k][a] = Inf
            end
        end
        sortperm!(ord[k], cst[k])
        for a in ord[k]
            cst[k][a] >= best[] && break
            u = levels[a]; U[k] = u
            dfs(X[k][a], k + 1, cst[k][a], u > 0 ? k : last, used + u * Δ)
        end
    end
    dfs(copy(x0), 1, 0.0, -10^6, 0.0)
    (U=copy(bestU), J=best[])
end

# ── closed loop ──────────────────────────────────────────────────────────────
# Measurement noise of the cell counts for bootstrap draw b: n×2 multiplicative factors (N, C), log-normal
# with mean 1 and coefficient of variation cv (control_config.toml [observation]). The seed depends on b
# only, so every policy simulated on draw b sees the same noise (paired comparisons).
function obs_noise(b, n; cv=Float64(CFG["observation"]["count_cv"]), seed=Int(CFG["observation"]["noise_seed"]))
    cv <= 0 && return ones(n, 2)
    σ = sqrt(log(1 + cv^2)); rng = Random.Xoshiro(seed + b)
    exp.(σ .* randn(rng, n, 2) .- σ^2 / 2)
end
# What the planner knows at decision k: its own model state xm (hidden states simulated under the applied
# inputs) with the cell counts replaced by the plant's counts × noise (exact counts when noise === nothing)
function observe!(xm, x, Pm, Pp, noise, k)
    xm[Pm.iN] = x[Pp.iN] * (noise === nothing ? 1.0 : noise[k, 1])
    xm[Pm.iC] = x[Pp.iC] * (noise === nothing ? 1.0 : noise[k, 2])
    xm
end

# Receding-horizon closed loop on the plant Pp, planning with Pm (the same model unless testing mismatch).
# The planner observes only N and C (observe!); the drug-effect chains, B and T are reconstructed by
# simulating Pm under the applied doses, as they are driven by the input alone. With Pm = Pp and no noise
# the reconstruction is exact, so the nominal closed loop equals replaying its own schedule (tested).
# Fixed look-ahead at every decision (shrink=true: the old horizon that ends with the run; comparison only).
# H: an Int (single intervals) or the look-ahead blocks (move blocking, discrete mode); continuous mode
# plans every interval of the same look-ahead.
function closed_loop(Pm::Plant, Pp::Plant, x0, Δ, n_dec, H, W::Weights; mode=:continuous, levels=Float64[],
                     shrink=false, cum=0.0, min_gap=0.0, u_init=0.5, noise=nothing)
    blk = mode == :continuous ? fill(1, sum(as_blocks(H))) : as_blocks(H); nb = length(blk)
    x = copy(x0); xm = copy(x0); applied = Float64[]; Uwarm = fill(min(u_init, W.umax), nb); used = 0.0; iters = Int[]
    for k in 1:n_dec
        observe!(xm, x, Pm, Pp, noise, k)
        bk = shrink ? cut_blocks(blk, n_dec - k + 1) : blk; Hk = length(bk)
        rem = cum > 0 ? max(cum - used, 0.0) : 0.0
        if mode == :continuous
            sol = sqp_solve(xm, Uwarm[1:Hk], Pm, Δ, W; cum=cum > 0 ? max(rem, 1e-9) : 0.0)
            push!(iters, sol.iters)
        else
            sol = discrete_solve(xm, Pm, Δ, bk, levels, W; min_gap=min_gap, cum=cum > 0 ? rem : 0.0, Uinit=Uwarm[1:Hk])
        end
        u = sol.U[1]; push!(applied, u); used += u * Δ
        x = F(x, u, Pp, Δ); xm = F(xm, u, Pm, Δ)
        Uwarm = vcat(sol.U[2:end], fill(sol.U[end], nb - Hk + 1))[1:nb]      # shift by one block (a bound only)
    end
    (U=applied, iters=iters, x=x)
end

# ── baseline policies (state feedback at decision times) ────────────────────
# noise: the same measurement noise as the closed-loop MPC (the policy sees N, C × noise[k, :])
function run_policy(policy, P::Plant, x0, Δ, n_dec; noise=nothing)
    x = copy(x0); U = Float64[]
    for k in 1:n_dec
        xo = noise === nothing ? x : (y = copy(x); y[P.iN] *= noise[k, 1]; y[P.iC] *= noise[k, 2]; y)
        u = policy(xo, k)
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

# Time to progression on a trajectory (progression: N + C ≥ B_prog):
#   first      the first day the burden reaches B_prog (NaN: never within the run)
#   loss       time to progression: the start of the first spell at or above B_prog that lasts at least
#              `sustain` days, or runs to the end of the run (NaN: none, i.e. not progressed). A brief
#              excursion that treatment brings back under the threshold within `sustain` days does not count
#              (treating hard at 80% confluence and backing off is fine), and a single dip under it after a
#              long spell above does not restart the clock.
#   days       days spent at or above B_prog
function time_to_progression(ts, tot, Bp; sustain=Float64(CFG["progression"]["sustain_days"]))
    first = findfirst(>=(Bp), tot)
    loss = NaN; i = 1; n = length(tot)
    while i <= n
        if tot[i] >= Bp
            j = i
            while j < n && tot[j+1] >= Bp
                j += 1
            end
            if j == n || ts[j] - ts[i] >= sustain
                loss = ts[i]; break
            end
            i = j + 1
        else
            i += 1
        end
    end
    days = sum((ts[i+1] - ts[i] for i in 1:length(ts)-1 if tot[i+1] >= Bp); init=0.0)
    (first=first === nothing ? NaN : ts[first], loss=loss, days=days)
end

# biologically interpretable outcomes of a schedule; days_over = days with N + C at or above B_cancer
function outcomes(P::Plant, x0, U, Δ, W::Weights)
    ts, xs = trajectory(x0, U, P, Δ)
    N = [x[P.iN] for x in xs]; C = [x[P.iC] for x in xs]; T = [x[P.iT] for x in xs]; tot = N .+ C
    over = findfirst(>=(W.Bc), tot)
    days_over = sum((ts[i+1] - ts[i] for i in 1:length(ts)-1 if tot[i+1] >= W.Bc); init=0.0)
    tp = time_to_progression(ts, tot, W.Bp)
    (N_final=N[end], C_final=C[end], total_final=tot[end], resistant_fraction=C[end] / max(tot[end], 1e-9),
     max_total=maximum(tot), max_tox=maximum(T), cum_dose=sum(U) * Δ, cancer_violation=over !== nothing,
     tox_violation=maximum(T) >= W.Btox, day_over_bound=over === nothing ? NaN : ts[over], days_over=days_over,
     ttp_first=tp.first, ttp=tp.loss, days_progressed=tp.days, ts=ts, N=N, C=C, T=T)
end
