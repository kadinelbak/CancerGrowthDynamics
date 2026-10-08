# Stage 2: universal differential equation (UDE). Every model is AUTONOMOUS (no clock-time terms),
# so it can be restarted from any state and used for optimal control (user rule, 2026-09-29).
#
#   dS/dt = r_S S (1 - (S + a_SR R)/K_S) - d_S(Q) S - kill_S(v, ρ) S
#   dR/dt = r_R R (1 - (R + a_RS S)/K_R) - d_R(Q) R - kill_R(v, ρ) R
#   dv/dt = γ (u - v)                      drug exposure (effect compartment driven by the dose u)
#   dQ/dt = -c Q (S + R) / 3500            resource left in low-resource wells (Q = 1 and constant in normal medium)
#   d_L(Q) = δ_L,low (1 - Q)               low-resource death rises as cells use up the resource
#   ρ = (S + R) / K_line                   relative density
# Compartments per line: cells + v + Q (≤ 3).
#
# variants
#   :none      no drug effect
#   :hill      kill = kmax_L v^n / (EC50_L^n + v^n)                      (mechanistic)
#   :symbolic  kill = closed form recovered from the NN (symbolic.jl)     (mechanistic)
#   :nn        kill = v softplus(NN_k(v, ρ, line))                        (UDE)
#   :nn_res    :nn, plus d_L(Q) = (1 - Q) softplus(NN_r(1 - Q, ρ, line))  (UDE, learned resource death)
#   :blackbox  per-capita growth = NN(ρ, v, line, low, Q)                (pure neural ODE)
# Competition a_SR, a_RS is fit separately on untreated co-culture (Step B); the pure neural ODE has no
# competition term and is fitted to the drug and untreated co-culture units together.
# Every well starts at t = 0 from its nominal seeding (estimated only where the first count is > 3^d-fold off).
# Everything is fixed-step RK4 in log-space-safe form and ForwardDiff-differentiable.

softplus(x) = log1p(exp(-abs(x))) + max(x, zero(x))

# ---------------- generic tanh MLP on a flat parameter vector ----------------
struct MLPSpec
    sizes::Vector{Int}
end
nparams(s::MLPSpec) = sum(s.sizes[i + 1] * (s.sizes[i] + 1) for i in 1:length(s.sizes) - 1)
function mlp_eval(s::MLPSpec, θ, x)
    k = 0; h = x
    for i in 1:length(s.sizes) - 1
        nin, nout = s.sizes[i], s.sizes[i + 1]
        W = reshape(view(θ, k + 1:k + nin * nout), nout, nin); k += nin * nout
        b = view(θ, k + 1:k + nout); k += nout
        z = W * h .+ b
        h = i < length(s.sizes) - 1 ? tanh.(z) : z
    end
    return h
end
function mlp_init(s::MLPSpec, rng; scale = 1.0)
    θ = Float64[]
    for i in 1:length(s.sizes) - 1
        nin, nout = s.sizes[i], s.sizes[i + 1]
        append!(θ, randn(rng, nout * nin) .* (scale / sqrt(nin)))
        append!(θ, zeros(nout))
    end
    return θ
end

const KILL_NET = MLPSpec([3, 8, 8, 1])
const BLACKBOX_NET = MLPSpec([5, 8, 8, 1])   # pure neural ODE comparator: [ρ, v, line, low, Q] (no time input)

"Generic allocation-free 8-8-1 tanh net with NIN inputs."
function net8(θ, off, x::SVector{NIN}) where {NIN}
    W1 = SMatrix{8, NIN}(ntuple(i -> θ[off + i], Val(8 * NIN)))
    b1 = SVector{8}(ntuple(i -> θ[off + 8 * NIN + i], Val(8)))
    h1 = tanh.(W1 * x + b1)
    o2 = off + 8 * NIN + 8
    W2 = SMatrix{8, 8}(ntuple(i -> θ[o2 + i], Val(64)))
    b2 = SVector{8}(ntuple(i -> θ[o2 + 64 + i], Val(8)))
    h2 = tanh.(W2 * h1 + b2)
    o3 = o2 + 72
    w3 = SVector{8}(ntuple(i -> θ[o3 + i], Val(8)))
    return dot(w3, h2) + θ[o3 + 9]
end

"Allocation-free evaluation of the 3-8-8-1 kill network (same layout as mlp_eval)."
function kill_net(θ, off, x1, x2, x3)
    W1 = SMatrix{8, 3}(ntuple(i -> θ[off + i], Val(24)))
    b1 = SVector{8}(ntuple(i -> θ[off + 24 + i], Val(8)))
    h1 = tanh.(W1 * SVector(x1, x2, x3) + b1)
    o2 = off + 32
    W2 = SMatrix{8, 8}(ntuple(i -> θ[o2 + i], Val(64)))
    b2 = SVector{8}(ntuple(i -> θ[o2 + 64 + i], Val(8)))
    h2 = tanh.(W2 * h1 + b2)
    o3 = o2 + 72
    w3 = SVector{8}(ntuple(i -> θ[o3 + i], Val(8)))
    return dot(w3, h2) + θ[o3 + 9]
end

# ---------------- parameter layout ----------------
# backbone (per line L in S,R): log r_L, log K_L; low-resource: log r_L,low, log δ_L,low; shared log c (resource use)
# drug: log γ (exposure filter); then either NN weights or Hill (log kmax_L, log EC50_L, log n)
const BACKBONE_NAMES = [:logr_S, :logK_S, :logr_R, :logK_R, :logr_S_low, :logδ_S_low,
                        :logr_R_low, :logδ_R_low, :logc, :logγ, :logfloor, :logfloor_low,
                        :log_r20_over_r30]      # growth-rate multiplier for 20k-seeded wells (density covariate)
const NB = length(BACKBONE_NAMES)
const HILL_NAMES = [:logkmax_S, :logEC50_S, :logkmax_R, :logEC50_R, :logn]

struct UDEModel
    kind::Symbol           # :nn, :hill, :none
    θ::Vector{Float64}     # backbone ++ kill parameters
    a::Vector{Float64}     # competition [a_SR, a_RS] (normal resource) ++ [a_SR, a_RS] (low resource)
end
nkill(kind) = kind == :nn ? nparams(KILL_NET) : kind == :nn_res ? 2 * nparams(KILL_NET) : kind == :hill ? length(HILL_NAMES) :
              kind == :blackbox ? nparams(BLACKBOX_NET) :
              kind == :symbolic ? length(SYMBOLIC_TERMS[].S) + length(SYMBOLIC_TERMS[].R) : 0

function init_theta(kind; rng = StableRNG(1))
    b = [log(0.9), log(3500.0), log(0.8), log(3500.0), log(0.9), log(0.3), log(0.8), log(0.5), log(0.3), log(0.24), log(5.0), log(10.0), 0.0]
    kind == :blackbox && return vcat(b, mlp_init(BLACKBOX_NET, rng; scale = 0.5))
    if kind == :nn_res
        k1 = mlp_init(KILL_NET, rng; scale = 0.5); k1[end] = -2.0
        k2 = mlp_init(KILL_NET, rng; scale = 0.5); k2[end] = -1.0
        return vcat(b, k1, k2)
    end
    k = kind == :nn ? mlp_init(KILL_NET, rng; scale = 0.5) :
        kind == :hill ? [log(0.5), log(1.2), log(0.1), log(3.0), log(3.0)] : Float64[]
    kind == :nn && (k[end] = -2.0)   # start with a small kill
    return vcat(b, k)
end

const KILL_MASK = Ref((1.0, 1.0, 1.0))

"Kill rate (1/day) for line (±1), exposure v and relative density ρ."
function kill_rate(kind, θ, v, ρ, line_sign)
    kind == :none && return zero(v * θ[1])
    if kind in (:nn, :nn_res)
        m = KILL_MASK[]     # input ablation (Stage 4); (1, 1, 1) = full model
        return v * softplus(kill_net(θ, NB, m[1] * v, m[2] * ρ, m[3] * line_sign))
    end
    kind == :symbolic && return symbolic_kill(θ, v, ρ, line_sign)
    h = view(θ, NB + 1:length(θ))
    kmax, ec50 = line_sign < 0 ? (exp(h[1]), exp(h[2])) : (exp(h[3]), exp(h[4]))
    n = exp(h[5]); vv = max(v, zero(v)) + 1e-9
    return kmax * vv^n / (ec50^n + vv^n)
end

"""
Right-hand side for one well: x = (S, R, v, Q). `t` only selects the dose u(t) from the schedule
(the control input); no rate depends on clock time.
"""
function well_rhs(kind, θ, a, low::Bool, sched::Schedule, t, x; d20::Bool = false)
    S, R, v = max(x[1], zero(x[1])), max(x[2], zero(x[2])), x[3]
    Q = clamp(x[4], zero(x[4]), one(x[4]))
    dQ = low ? -exp(θ[9]) * Q * (S + R) / 3500 : zero(exp(θ[9]) * Q * S)
    if kind == :blackbox
        ρ = (S + R) / 3500
        lw = low ? 1.0 : -1.0
        gS = net8(θ, NB, SVector(ρ, v, -1.0, lw, Q))
        gR = net8(θ, NB, SVector(ρ, v, 1.0, lw, Q))
        return (gS * S, gR * R, exp(θ[10]) * (dose_at(sched, t) - v), dQ)
    end
    rS, KS, rR, KR = exp(θ[1]), exp(θ[2]), exp(θ[3]), exp(θ[4])
    γ = exp(θ[10])
    δS = δR = zero(rS * Q)
    if low
        rS, rR = exp(θ[5]), exp(θ[7])
        used = 1 - Q
        if kind == :nn_res
            off = NB + nparams(KILL_NET)
            δS = used * softplus(kill_net(θ, off, used, (S + R) / KS, -1.0))
            δR = used * softplus(kill_net(θ, off, used, (S + R) / KR, 1.0))
        else
            δS, δR = exp(θ[6]) * used, exp(θ[8]) * used
        end
    end
    if d20                                  # density covariate: 20k wells grow at r × exp(θ[13])
        m20 = exp(θ[13]); rS *= m20; rR *= m20
    end
    kS = kill_rate(kind, θ, v, (S + R) / KS, -1.0)
    kR = kill_rate(kind, θ, v, (S + R) / KR, 1.0)
    aSR, aRS = low ? (a[3], a[4]) : (a[1], a[2])
    dS = rS * S * (1 - (S + aSR * R) / KS) - δS * S - kS * S
    dR = rR * R * (1 - (R + aRS * S) / KR) - δR * R - kR * R
    dv = γ * (dose_at(sched, t) - v)
    return (dS, dR, dv, dQ)
end

"""
Integrate a well from t0 with state (S0, R0, v0, Q0) and return (S, R, v, Q) at `tobs` (sorted, ≥ t0).
Fixed-step RK4; steps are aligned so dose switches at integer days are hit exactly.
"""
function simulate_well(kind, θ, a, low, sched, t0, S0, R0, v0, tobs; dt = 0.1, d20::Bool = false, Q0 = 1.0)
    T = promote_type(eltype(θ), typeof(S0), typeof(R0), typeof(v0), typeof(Q0), eltype(a), eltype(sched.doses))
    x = (T(S0), T(R0), T(v0), T(Q0))
    t = Float64(t0)
    out = Vector{NTuple{4, T}}(undef, length(tobs))
    j = 1
    while j <= length(tobs) && tobs[j] <= t + 1e-9
        out[j] = x; j += 1
    end
    f(tt, xx) = well_rhs(kind, θ, a, low, sched, tt, xx; d20)
    while j <= length(tobs)
        h = min(dt, tobs[j] - t)
        k1 = f(t, x)
        k2 = f(t + h / 2, x .+ (h / 2) .* k1)
        k3 = f(t + h / 2, x .+ (h / 2) .* k2)
        k4 = f(t + h, x .+ h .* k3)
        x = x .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        x = (max(x[1], T(1e-6)), max(x[2], T(1e-6)), x[3], clamp(x[4], zero(T), one(T)))
        t += h
        while j <= length(tobs) && tobs[j] <= t + 1e-9
            out[j] = x; j += 1
        end
    end
    return out
end

# ---------------- data -> fitting units ----------------
const SIGMA_FLOOR = 0.15   # log-scale noise floor (plate-level artefacts, e.g. shared day-8 dips)

"""
A fitting unit is one well: its members' observations (log scale) on the measured days (≥ 1). Every well
starts at t = 0 from its nominal seeding (day 0 is not a measurement); a member whose first count is more
than 3^d-fold off the nominal value has an estimated starting value (`free`), bounded as in zoo_data.jl.
"""
struct Unit
    members::Vector{Trajectory}
    low::Bool
    sched::Schedule
    tobs::Vector{Float64}
    y::Vector{Vector{Float64}}      # per member, log obs on tobs (NaN = missing)
    w::Vector{Vector{Float64}}      # 1/σ²
    l0nom::Vector{Float64}          # per member log nominal starting value
    free::Vector{Bool}              # per member: starting value estimated
    lo::Vector{Float64}             # log bounds of an estimated starting value
    hi::Vector{Float64}
    z0::Vector{Float64}             # initial unconstrained values of the free members (start = √(n0 · first count))
end

function make_unit(members)
    sched = members[1].schedule
    tobs = sort(unique(vcat((tr.days for tr in members)...)))
    y = Vector{Float64}[]; w = Vector{Float64}[]
    for tr in members
        d = Dict(zip(tr.days, eachindex(tr.days)))
        push!(y, [haskey(d, t) ? log(max(tr.mean[d[t]], 1e-1)) : NaN for t in tobs])
        σ = [haskey(d, t) ? sqrt(SIGMA_FLOOR^2 + (tr.sem[d[t]] / max(tr.mean[d[t]], 1e-1))^2) : Inf for t in tobs]
        push!(w, 1 ./ σ .^ 2)
    end
    l0nom = [log(tr.n0) for tr in members]
    free = [tr.n0_free for tr in members]
    lo = [log(n0_bounds(tr)[1]) for tr in members]; hi = [log(n0_bounds(tr)[2]) for tr in members]
    z0 = Float64[]
    for (i, tr) in enumerate(members)
        free[i] || continue
        l = clamp((l0nom[i] + log(tr.mean[1])) / 2, lo[i] + 0.01 * (hi[i] - lo[i]), hi[i] - 0.01 * (hi[i] - lo[i]))
        q = (l - lo[i]) / (hi[i] - lo[i])
        push!(z0, log(q / (1 - q)))
    end
    return Unit(members, members[1].resource == :low, sched, tobs, y, w, l0nom, free, lo, hi, z0)
end

make_units(trajs) = [make_unit(w) for w in wells(trajs)]
is20k(u) = u.members[1].seeding == 20
nfree(u::Unit) = count(u.free)

"Per-member log starting values: nominal, or the bounded transform of `z` for the estimated members."
function unit_l0(u::Unit, z)
    T = promote_type(Float64, eltype(z))
    l = Vector{T}(undef, length(u.members)); k = 0
    for i in eachindex(u.members)
        if u.free[i]
            k += 1
            l[i] = u.lo[i] + (u.hi[i] - u.lo[i]) / (1 + exp(-z[k]))
        else
            l[i] = u.l0nom[i]
        end
    end
    return l
end

"Predicted log values per member on `u.tobs` from log starting values `l0` at t = 0 (v = 0, Q = 1)."
function predict_unit(kind, θ, a, u::Unit, l0)
    S0 = zero(eltype(l0)); R0 = zero(eltype(l0))
    for (i, tr) in enumerate(u.members)
        tr.line == :S ? (S0 = exp(l0[i])) : (R0 = exp(l0[i]))
    end
    xs = simulate_well(kind, θ, a, u.low, u.sched, 0.0, S0, R0, 0.0, u.tobs; d20 = is20k(u))
    floor = exp(u.low ? θ[12] : θ[11])   # detection floor (image-analysis background), count scale; separate for low resource
    return [[log((tr.line == :S ? x[1] : x[2]) + floor) for x in xs] for tr in u.members]
end

"Weighted SSE (log scale) for a unit."
function unit_sse(kind, θ, a, u::Unit, l0)
    pred = predict_unit(kind, θ, a, u, l0)
    s = zero(eltype(θ)) + zero(eltype(l0)) + zero(eltype(a))
    for i in eachindex(u.members), k in eachindex(u.tobs)
        isnan(u.y[i][k]) && continue
        s += u.w[i][k] * (u.y[i][k] - pred[i][k])^2
    end
    return s
end

# ---------------- Step A: fit backbone + kill on :drug units ----------------
function _unpack(p, nθ, units)
    θ = view(p, 1:nθ); k = nθ; zs = []
    for u in units
        n = nfree(u); push!(zs, view(p, k + 1:k + n)); k += n
    end
    return θ, zs
end

"Weak priors that keep backbone rates in a biological range (identical for every model kind)."
function backbone_prior(θ)
    return sum(abs2, (θ[[1, 3, 5, 7]] .- log(0.8)) ./ 1.5) + sum(abs2, (θ[[2, 4]] .- log(3500.0)) ./ 1.0) +
           abs2((θ[10] - log(0.24)) / 1.0) + abs2((θ[9] - log(0.3)) / 1.5) + abs2(θ[11] - log(5.0)) + abs2(θ[12] - log(10.0)) + abs2(θ[13])
end

"""
fit_drug_model(kind, units; ...) jointly fits backbone + drug term + the estimated starting values.
Adam warm-up then L-BFGS, gradients via ForwardDiff. Returns (θ, l0s, loss).
"""
function fit_drug_model(kind, units; θinit = init_theta(kind), l2 = 1e-3, adam_steps = 600, lbfgs_iters = 400,
                        lr = 0.02, verbose = false, weights = ones(length(units)))
    nθ = length(θinit)
    p0 = vcat(θinit, reduce(vcat, [u.z0 for u in units]; init = Float64[]))
    a0 = ones(4)
    function loss(p)
        θ, zs = _unpack(p, nθ, units)
        s = zero(eltype(p))
        for (u, z, wu) in zip(units, zs, weights)
            s += wu * unit_sse(kind, θ, a0, u, unit_l0(u, z))
        end
        kind in (:nn, :nn_res, :blackbox) && (s += l2 * sum(abs2, view(θ, NB + 1:nθ)) * length(units))
        return s + backbone_prior(θ)
    end
    cfg = ForwardDiff.GradientConfig(loss, p0, ForwardDiff.Chunk{min(12, length(p0))}())
    g!(G, p) = ForwardDiff.gradient!(G, loss, p, cfg)
    p = copy(p0); m = zero(p); v = zero(p); G = zero(p)
    for k in 1:adam_steps
        g!(G, p)
        all(isfinite, G) || (G .= 0)
        @. m = 0.9m + 0.1G; @. v = 0.999v + 0.001G^2
        @. p -= lr * (m / (1 - 0.9^k)) / (sqrt(v / (1 - 0.999^k)) + 1e-8)
        verbose && k % 100 == 0 && @info "adam" k loss(p)
    end
    res = Optim.optimize(loss, g!, p, LBFGS(), Optim.Options(iterations = lbfgs_iters, g_tol = 1e-6))
    p = Optim.minimizer(res)
    θ, zs = _unpack(p, nθ, units)
    return (θ = collect(θ), l0s = [unit_l0(u, collect(z)) for (u, z) in zip(units, zs)], loss = Optim.minimum(res))
end

# ---------------- Step B: competition on untreated co-culture ----------------
"""
Fit competition a (normal and low resource) and re-estimate the low-resource detection floor on
untreated co-culture; drug/backbone parameters stay fixed. Returns the calibrated UDEModel.
"""
function calibrate(kind, θ, units)
    function loss(p)
        a = exp.(p[1:4])
        θc = vcat(θ[1:11], p[5], θ[13:end])
        s = zero(eltype(p)); k = 5
        for u in units
            n = nfree(u); z = p[k + 1:k + n]; k += n
            s += unit_sse(kind, θc, a, u, unit_l0(u, z))
        end
        return s + sum(abs2, p[1:4] ./ 2) + abs2(p[5] - θ[12])   # weak priors around a = 1 and the drug-stage floor
    end
    p0 = vcat(zeros(4), θ[12], reduce(vcat, [u.z0 for u in units]; init = Float64[]))
    g!(G, p) = ForwardDiff.gradient!(G, loss, p)
    res = Optim.optimize(loss, g!, p0, LBFGS(), Optim.Options(iterations = 300))
    p = Optim.minimizer(res)
    θc = copy(θ); θc[12] = p[5]
    return UDEModel(kind, θc, exp.(p[1:4]))
end

# ---------------- evaluation (forward prediction from t = 0, no refit of the model) ----------------
const START_CACHE = Dict{Tuple{UInt, String}, Vector{Float64}}()
const START_LOCK = ReentrantLock()
"""
Log starting values of a TRAINING unit under a fitted model: the nominal seeding, and for members whose
starting value is estimated (first count > 3^d-fold off nominal), the value that best fits that member's
own data with every model parameter frozen. Held-out wells use the nominal seeding instead (`nominal = true`).
"""
function start_values(model::UDEModel, u::Unit)
    nfree(u) == 0 && return copy(u.l0nom)
    key = (hash((model.kind, model.θ, model.a, u.tobs, u.y)), u.members[1].group)
    hit = lock(() -> get(START_CACHE, key, nothing), START_LOCK)
    hit === nothing || return copy(hit)
    l0 = _start_values(model, u)
    lock(() -> (START_CACHE[key] = l0), START_LOCK)
    return copy(l0)
end
function _start_values(model::UDEModel, u::Unit)
    f(z) = unit_sse(model.kind, model.θ, model.a, u, unit_l0(u, z))
    g!(G, z) = ForwardDiff.gradient!(G, f, z)
    best = nothing
    for z0 in (u.z0, fill(-1.0, nfree(u)), fill(1.0, nfree(u)))
        r = Optim.optimize(f, g!, copy(z0), LBFGS(), Optim.Options(iterations = 200))
        (best === nothing || Optim.minimum(r) < Optim.minimum(best)) && (best = r)
    end
    return unit_l0(u, Optim.minimizer(best))
end

"""
Per-point log observations and predictions of a unit, forward from t = 0 (no refit). `nominal = true`
(every held-out well) starts every member from its nominal seeding, as the zoo's validation does
(02_mechanical_automatic_package/validation.toml): held-out data are never used to set a starting value.
"""
function forward_residuals(model::UDEModel, u::Unit; nominal::Bool = false, l0 = nominal ? copy(u.l0nom) : start_values(model, u))
    pred = predict_unit(model.kind, model.θ, model.a, u, l0)
    rows = DataFrame(traj_id = String[], day = Float64[], obs = Float64[], pred = Float64[])
    for (i, tr) in enumerate(u.members), k in findall(!isnan, u.y[i])
        push!(rows, (tr.id, u.tobs[k], u.y[i][k], pred[i][k]))
    end
    return rows
end

"""
Rolling one-step (next-day) prediction from each observed day (UDE analogue of Stage 1). Cells restart
from the observation; the unobserved states v and Q come from the model's own forward run.
"""
function nextday_residuals(model::UDEModel, u::Unit; nominal::Bool = false)
    out = Float64[]
    l0 = nominal ? copy(u.l0nom) : start_values(model, u)
    S0f = R0f = 0.0
    for (i, tr) in enumerate(u.members)
        tr.line == :S ? (S0f = exp(l0[i])) : (R0f = exp(l0[i]))
    end
    xs = simulate_well(model.kind, model.θ, model.a, u.low, u.sched, 0.0, S0f, R0f, 0.0, u.tobs; d20 = is20k(u))
    for (k, t) in enumerate(u.tobs)
        k == length(u.tobs) && break
        all(!isnan(u.y[i][k]) && !isnan(u.y[i][k + 1]) for i in eachindex(u.members)) || continue
        u.tobs[k + 1] - t == 1 || continue
        S0 = R0 = 0.0
        for (i, tr) in enumerate(u.members)
            val = exp(u.y[i][k])
            tr.line == :S ? (S0 = val) : (R0 = val)
        end
        x = simulate_well(model.kind, model.θ, model.a, u.low, u.sched, t, S0, R0, xs[k][3], [t + 1]; d20 = is20k(u), Q0 = xs[k][4])[1]
        for (i, tr) in enumerate(u.members)
            p = log((tr.line == :S ? x[1] : x[2]) + exp(model.θ[u.low ? 12 : 11]))
            push!(out, p - u.y[i][k + 1])
        end
    end
    return out
end

# ---------------- ensembles, evaluation helpers, persistence ----------------
"""
Random-seed replicates: the same data and unit weights, `n` different random initialisations of the
network weights (seed i uses StableRNG(1000 + i)). Used for every neural model; the spread over seeds is
the optimisation / initialisation uncertainty.
"""
function fit_seeds(kind, units; n = 5, θwarm = nothing, kw...)
    members = Vector{Any}(undef, n)
    Threads.@threads for i in 1:n
        θ0 = init_theta(kind; rng = StableRNG(1000 + i))
        θwarm === nothing || (θ0[1:NB] .= θwarm[1:NB])
        members[i] = fit_drug_model(kind, units; θinit = θ0, kw...)
    end
    return members
end

"Bayesian-bootstrap ensemble: member 1 uses unit weights 1, others Dirichlet(1) weights × n."
function fit_ensemble(kind, units; n = 8, θwarm = nothing, kw...)
    members = Vector{Any}(undef, n)
    Threads.@threads for i in 1:n
        rng = StableRNG(1000 + i)
        w = i == 1 ? ones(length(units)) : (e = -log.(rand(rng, length(units))); e ./ sum(e) .* length(units))
        θ0 = init_theta(kind; rng)
        θwarm === nothing || (θ0[1:NB] .= θwarm[1:NB]; kind == :hill && (θ0 .= θwarm))
        members[i] = fit_drug_model(kind, units; θinit = θ0, weights = w, kw...)
    end
    return members
end

"Mean log prediction across an ensemble of UDEModels, plus 10/90% band."
function ensemble_forward(models::Vector{UDEModel}, u::Unit; nominal::Bool = false)
    rs = [forward_residuals(m, u; nominal) for m in models]
    base = rs[1][:, [:traj_id, :day, :obs]]
    P = reduce(hcat, [r.pred for r in rs])
    base.pred = vec(mean(P; dims = 2))
    base.lo = [quantile(P[i, :], 0.1) for i in 1:size(P, 1)]
    base.hi = [quantile(P[i, :], 0.9) for i in 1:size(P, 1)]
    return base
end

rmse(x) = sqrt(mean(abs2, x))

function save_models(path, dict)
    open(path, "w") do io
        JSON3.write(io, Dict(k => [Dict("kind" => String(m.kind), "theta" => m.θ, "a" => m.a) for m in v] for (k, v) in dict))
    end
end
function load_models(path)
    raw = JSON3.read(read(path, String))
    return Dict(String(k) => [UDEModel(Symbol(m.kind), Float64.(m.theta), Float64.(m.a)) for m in v] for (k, v) in raw)
end

"""
Dynamic parallel map: every thread pulls the next task from a shared queue as soon as it is free (unlike
`Threads.@threads`, which splits the list into fixed chunks), so long and short fits balance out.
"""
function pmap_dynamic(f, xs)
    out = Vector{Any}(undef, length(xs))
    ch = Channel{Int}(length(xs))
    foreach(i -> put!(ch, i), eachindex(xs)); close(ch)
    @sync for _ in 1:Threads.nthreads()
        Threads.@spawn for i in ch
            out[i] = f(xs[i])
        end
    end
    return out
end
