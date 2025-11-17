
# 4) Utilities and models
struct FitResult
    r::Float64
    K::Float64
    SSR::Float64
    BIC::Float64
end

function logistic!(du, u, p, t)
    r, K = p
    du[1] = r * u[1] * (1 - u[1] / K)
end

function find_col(df::DataFrame; must_include::Vector{String}=String[])
    for c in names(df)
        s = lowercase(String(c))
        if all(substr -> occursin(substr, s), must_include)
            return Symbol(c)
        end
    end
    return nothing
end

function load_day_averages(path::AbstractString)
    df = CSV.read(path, DataFrame)
    daycol = (:Day in names(df)) ? :Day : find_col(df; must_include=["day"])
    daycol === nothing && error("No day column in $(basename(path))")

    val_candidates = [
        Symbol("Mean Cells"), Symbol("Mean_Cells"), Symbol("mean_cells"), Symbol("day_mean_cells"),
        Symbol("Mean_Cell_Count"), Symbol("MeanCells"), Symbol("Mean cells")
    ]
    valcol = nothing
    for cand in val_candidates
        if cand in names(df)
            valcol = cand
            break
        end
    end
    if valcol === nothing
        valcol = find_col(df; must_include=["mean", "cell"])
    end
    valcol === nothing && error("No mean cell column in $(basename(path))")

    errcol = find_col(df; must_include=["sem", "cell"])
    if isnothing(errcol)
        errcol = find_col(df; must_include=["sd", "cell"])
    end

    x = Float64.(df[!, daycol])
    y = Float64.(df[!, valcol])
    yerr = isnothing(errcol) ? nothing : Float64.(df[!, errcol])

    perm = sortperm(x)
    x = x[perm]; y = y[perm]; yerr = isnothing(yerr) ? nothing : yerr[perm]
    return x, y, yerr
end

function infer_density(path::AbstractString)
    for part in reverse(splitpath(String(path)))
        lp = lowercase(part)
        if lp == "20k" || lp == "30k"
            return lp
        end
    end
    return ""
end

function parse_coculture_meta(path::AbstractString)
    name = basename(path)
    mix_match = match(r"measure_([0-9]+-[0-9]+)", name)
    ratio = mix_match === nothing ? "" : (isempty(mix_match.captures) ? "" : mix_match.captures[1])
    variant = occursin("cis", lowercase(name)) ? "A2780cis" : "A2780"
    density = infer_density(path)
    return (ratio=ratio, variant=variant, density=density)
end

function fit_logistic(x::Vector{<:Real}, y::Vector{<:Real}; r_bounds::Tuple{Float64,Float64}=(0.0, 2.0), K_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing, density::String="")
    x = collect(Float64.(x)); y = collect(Float64.(y))
    day0_cells = density == "30k" ? 100.0 : density == "20k" ? 67.0 : max(y[1], eps())

    y_max = maximum(y)
    y_norm = y ./ y_max

    tspan = (x[1], x[end])
    u0 = [day0_cells / y_max]

    k_hi = isnothing(K_bounds) ? 5.0 : K_bounds[2] / y_max
    k_lo = isnothing(K_bounds) ? 0.1 : K_bounds[1] / y_max
    bounds = [r_bounds, (k_lo, k_hi)]
    solver = Rosenbrock23()

    best_loss = Inf
    best_p = nothing
    best_sol = nothing
    for _ in 1:10
        p0 = [rand(Uniform(r_bounds...)), rand(Uniform(bounds[2][1], bounds[2][2]))]
        prob = ODEProblem(logistic!, u0, tspan, p0)
        obj = build_loss_objective(prob, solver, L2Loss(x, y_norm), Optimization.AutoForwardDiff())
        res = bboptimize(obj; SearchRange=bounds, MaxTime=15.0, TraceMode=:silent)
        pI = best_candidate(res)
        solI = solve(remake(prob, p=pI), solver; saveat=x, reltol=1e-9, abstol=1e-9)
        pred = getindex.(solI.u, 1)
        loss = sum(abs2.(y_norm .- pred))
        if loss < best_loss
            best_loss = loss
            best_p = pI
            best_sol = solI
        end
    end

    @assert best_p !== nothing "Fit failed to converge"
    r, k_norm = best_p
    K = k_norm * y_max

    sol_original = solve(ODEProblem(logistic!, [day0_cells], tspan, [r, K]), solver; saveat=x)
    pred_original = getindex.(sol_original.u, 1)
    ssr = sum(abs2.(y .- pred_original))
    n = length(y)
    bic = n * log(ssr / n) + 2 * log(n)

    return FitResult(r, K, ssr, bic), sol_original, (r=r, K=K, SSR=ssr, BIC=bic)
end

function plot_data_and_fit(x, y, sol; title_str::String="Logistic Fit", yerr=nothing)
    if yerr === nothing
        p = scatter(x, y; label="Data", xlabel="Day", ylabel="Cells", title=title_str, markersize=4)
    else
        p = scatter(x, y; yerr=yerr, label="Data", xlabel="Day", ylabel="Cells", title=title_str, markersize=4)
    end
    plot!(p, sol.t, getindex.(sol.u, 1); label="Logistic Model", linewidth=2, color=:red)
    display(p)
    return p
end

# Temporary test on synthetic data
_test_x = [0.0, 1.0, 2.0, 3.0]
_test_y = [50.0, 80.0, 140.0, 230.0]
_test_res, _test_sol, _ = fit_logistic(_test_x, _test_y; density="20k")
@assert _test_res.r > 0 && _test_res.K > maximum(_test_y)

