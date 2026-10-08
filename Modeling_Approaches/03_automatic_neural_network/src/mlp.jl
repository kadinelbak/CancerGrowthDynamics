# Small one-hidden-layer MLP with hand-written backprop + Adam (Stage 1 next-day model).

struct MLP1
    W1::Matrix{Float64}; b1::Vector{Float64}; w2::Vector{Float64}; b2::Float64
    mu::Vector{Float64}; sd::Vector{Float64}          # feature standardisation
end

function (m::MLP1)(X::AbstractMatrix)                 # X: features x samples (raw scale)
    Z = (X .- m.mu) ./ m.sd
    H = tanh.(m.W1 * Z .+ m.b1)
    return vec(m.w2' * H) .+ m.b2
end

function fit_mlp1(X, y; hidden = 12, lambda = 1e-3, steps = 3000, lr = 0.01, rng = StableRNG(1))
    p, n = size(X)
    mu = vec(mean(X; dims = 2)); sd = vec(std(X; dims = 2)); sd[sd .< 1e-8] .= 1.0
    Z = (X .- mu) ./ sd
    W1 = randn(rng, hidden, p) ./ sqrt(p); b1 = zeros(hidden); w2 = randn(rng, hidden) ./ sqrt(hidden); b2 = [mean(y)]
    params = (W1, b1, w2, b2)
    m1 = map(zero, params); m2 = map(zero, params)
    β1, β2, ϵ = 0.9, 0.999, 1e-8
    for k in 1:steps
        H = tanh.(W1 * Z .+ b1)
        e = (vec(w2' * H) .+ b2[1] .- y) .* (2 / n)
        dH = (w2 * e') .* (1 .- H .^ 2)
        g = (dH * Z' .+ 2lambda .* W1, vec(sum(dH; dims = 2)), H * e .+ 2lambda .* w2, [sum(e)])
        for (P, G, M, V) in zip(params, g, m1, m2)
            @. M = β1 * M + (1 - β1) * G
            @. V = β2 * V + (1 - β2) * G^2
            @. P -= lr * (M / (1 - β1^k)) / (sqrt(V / (1 - β2^k)) + ϵ)
        end
    end
    return MLP1(W1, b1, w2, b2[1], mu, sd)
end

struct MLPEnsemble
    members::Vector{MLP1}
end
(e::MLPEnsemble)(X) = mean(m(X) for m in e.members)
fit_mlp_ensemble(X, y; n = 5, seed0 = 100, kw...) = MLPEnsemble([fit_mlp1(X, y; rng = StableRNG(seed0 + i), kw...) for i in 1:n])
