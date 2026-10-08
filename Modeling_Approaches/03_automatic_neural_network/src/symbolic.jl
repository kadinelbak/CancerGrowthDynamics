# Stage 3: turn the learned NN kill term back into an equation.
# 1. sample the UDE ensemble's kill(v, ρ) per line on the region the data actually visit;
# 2. sparse regression over a small library (exhaustive best subset ≤ 3 terms, plus STLSQ sweep);
# 3. best form per number of terms (the front);
# 4. refit (A2780 form, A2780cis form) candidates as ODEs on the data, keep the lowest BIC (run_04_symbolic.jl).

"Library term: name + function of (v, ρ)."
struct Term
    name::String
    f::Function
end

hill(v, c, n) = (vv = max(v, zero(v)) + 1e-9; vv^n / (c^n + vv^n))

const LIBRARY = vcat(
    [Term("v", (v, ρ) -> v), Term("v²", (v, ρ) -> v^2), Term("v³", (v, ρ) -> v^3),
     Term("v·ρ", (v, ρ) -> v * ρ), Term("v²·ρ", (v, ρ) -> v^2 * ρ), Term("v³·ρ", (v, ρ) -> v^3 * ρ)],
    [Term("H(v; $(c), $(n))", (v, ρ) -> hill(v, c, n)) for c in (0.8, 1.0, 1.2, 1.5) for n in (2, 4)],
)

"States (v, ρ) visited by the fitted model along the training drug units, per line."
function visited_states(model::UDEModel, units)
    pts = Dict(-1.0 => NTuple{2, Float64}[], 1.0 => NTuple{2, Float64}[])
    for u in units
        l0 = start_values(model, u)
        S0 = R0 = 0.0
        for (i, tr) in enumerate(u.members)
            tr.line == :S ? (S0 = exp(l0[i])) : (R0 = exp(l0[i]))
        end
        ts = collect(0.0:0.25:maximum(u.tobs))
        xs = simulate_well(model.kind, model.θ, model.a, u.low, u.sched, 0.0, S0, R0, 0.0, ts; d20 = is20k(u))
        KS, KR = exp(model.θ[2]), exp(model.θ[4])
        for x in xs
            x[3] > 1e-3 || continue
            u.members[1].line == :S && push!(pts[-1.0], (x[3], (x[1] + x[2]) / KS))
            u.members[1].line == :R && push!(pts[1.0], (x[3], (x[1] + x[2]) / KR))
        end
    end
    return pts
end

"Grid over the visited box (v up to the max visited exposure) and ensemble-mean NN kill."
function kill_samples(models::Vector{UDEModel}, units; nv = 25, nρ = 12)
    out = Dict{Float64, NamedTuple}()
    for line in (-1.0, 1.0)
        pts = reduce(vcat, [visited_states(m, units)[line] for m in models[1:min(3, end)]])
        vmax = maximum(first.(pts)); ρlo, ρhi = extrema(last.(pts))
        V = Float64[]; P = Float64[]; K = Float64[]; Ksd = Float64[]
        for v in range(0.0, vmax; length = nv), ρ in range(ρlo, ρhi; length = nρ)
            ks = [kill_rate(:nn, m.θ, v, ρ, line) for m in models]
            push!(V, v); push!(P, ρ); push!(K, mean(ks)); push!(Ksd, std(ks))
        end
        out[line] = (v = V, ρ = P, k = K, ksd = Ksd, vmax = vmax, ρrange = (ρlo, ρhi))
    end
    return out
end

design(terms, v, ρ) = reduce(hcat, [[t.f(vi, ρi) for (vi, ρi) in zip(v, ρ)] for t in terms])

"Exhaustive best subset up to `maxterms`; returns DataFrame sorted by (nterms, rmse)."
function best_subsets(samples; maxterms = 3)
    rows = DataFrame(nterms = Int[], terms = Vector{Int}[], formula = String[], coefs = Vector{Float64}[], rmse = Float64[], r2 = Float64[])
    Φall = design(LIBRARY, samples.v, samples.ρ)
    y = samples.k
    sst = sum(abs2, y .- mean(y))
    function subsets(n, k)
        k == 0 && return [Int[]]
        k > n && return Vector{Int}[]
        return vcat([vcat(s, n) for s in subsets(n - 1, k - 1)], subsets(n - 1, k))
    end
    push!(rows, (0, Int[], "0", Float64[], sqrt(mean(abs2, y)), 1 - sum(abs2, y) / sst))
    for k in 1:maxterms, idx in subsets(length(LIBRARY), k)
        idx = sort(idx)
        Φ = Φall[:, idx]
        c = Φ \ y
        r = y .- Φ * c
        push!(rows, (k, idx, join(["$(round(ci; sigdigits = 3))·$(LIBRARY[j].name)" for (ci, j) in zip(c, idx)], " + "),
                     c, sqrt(mean(abs2, r)), 1 - sum(abs2, r) / sst))
    end
    sort!(rows, [:nterms, :rmse])
    return rows
end

"Sequentially thresholded least squares (SINDy) on normalised columns, swept over λ."
function stlsq_sweep(samples; lambdas = exp10.(range(-3, 0; length = 25)))
    Φ = design(LIBRARY, samples.v, samples.ρ)
    sc = [maximum(abs, Φ[:, j]) for j in 1:size(Φ, 2)]; sc[sc .== 0] .= 1
    Z = Φ ./ sc'
    y = samples.k
    rows = DataFrame(lambda = Float64[], nterms = Int[], formula = String[], rmse = Float64[])
    for λ in lambdas
        active = trues(size(Z, 2))
        c = zeros(size(Z, 2))
        for _ in 1:20
            any(active) || break
            c .= 0
            c[active] = Z[:, active] \ y
            small = abs.(c) .< λ * maximum(abs, y)
            newactive = active .& .!small
            newactive == active && break
            active = newactive
        end
        idx = findall(active)
        coefs = c[idx] ./ sc[idx]
        push!(rows, (λ, length(idx), join(["$(round(ci; sigdigits = 3))·$(LIBRARY[j].name)" for (ci, j) in zip(coefs, idx)], " + "),
                     sqrt(mean(abs2, y .- Φ[:, idx] * coefs))))
    end
    return unique(rows, :formula)
end

# ---------------- closed form as an ODE (kind :symbolic) ----------------
# θ layout: backbone (NB) ++ coefficients for S terms ++ coefficients for R terms.
const SYMBOLIC_TERMS = Ref((S = Int[], R = Int[]))

function symbolic_kill(θ, v, ρ, line_sign)
    terms = line_sign < 0 ? SYMBOLIC_TERMS[].S : SYMBOLIC_TERMS[].R
    off = NB + (line_sign < 0 ? 0 : length(SYMBOLIC_TERMS[].S))
    k = zero(v * θ[1])
    for (i, j) in enumerate(terms)
        k += θ[off + i] * LIBRARY[j].f(v, ρ)
    end
    return max(k, zero(k))
end
