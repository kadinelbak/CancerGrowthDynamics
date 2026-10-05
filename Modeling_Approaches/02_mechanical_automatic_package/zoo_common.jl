# Data loading and fitting helpers for the autonomous model zoo (run_zoo_fits.jl).
# Standalone: it does not include run_new_data_staged_fits.jl. The data come from the tidy table
# (zoo_data.jl), which reproduces that script's files, excluded points and day-0 shifts exactly.

import Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "stage_names.jl"))   # thesis stage names for every page (1, 2a, 2b, 3a, 3b, 4)
using GrowthParameterEstimation
using CSV, DataFrames, XLSX, Dates, Printf, TOML, SHA
using OrdinaryDiffEq: ODEProblem, solve, Tsit5
using LinearAlgebra: pinv, diag

const ROOT     = normpath(joinpath(@__DIR__, "..", ".."))
const EC50     = 1.0                        # µM, fixed
const FBS_LR   = 0.5
const U_R      = (10 - FBS_LR) / (10 - 0.05)   # low-resource level: 0 in 10% FBS, 0.955 in 0.5% FBS

# ── data ───────────────────────────────────────────────────────────────────────
# Every data set comes from the tidy table data/measurements.csv (built from data/sources.toml by
# data/build_tidy.jl), selected by its tags in zoo_data.jl: files, sheets, doses, time origins and
# QC exclusions are configuration, not code.
include(joinpath(@__DIR__, "zoo_data.jl"))
const COCULTURE_DOSE = maximum(e.dose for e in COCULTURE_LAYOUT)   # µM in the treated (Ce1) cocultures

const XMAX = 14.0
const YMAX = ceil(maximum(vcat([u.y for u in values(untreated)]..., [y for t in values(treated) for (_, y) in t.data]...,
                               [d.S[2] for d in co_data]..., [d.R[2] for d in co_data]...)) * 1.10; sigdigits=2)

# ── fitting helpers (used by bootstrap_tools.jl too) ────────────────────────────

function dense(model, u0, p, tmax; dt=0.1)
    ts = collect(0.0:dt:tmax)
    sol = solve(ODEProblem(model, u0, (0.0, tmax), p), Tsit5(); saveat=ts, reltol=1e-8, abstol=1e-8)
    ts, sol
end

# per-point residual scales of a dataset (1 when unweighted, see point_scale in zoo_data.jl)
dweights(d) = hasproperty(d, :w) ? d.w : Iterators.repeated(1.0, length(d.x))
# datasets as GrowthParameterEstimation takes them: a weighted dataset becomes one dataset per point,
# each with its residual_scale (the package solves the model once per evaluation, so this is cheap)
gpe_datasets(ds) = any(d -> hasproperty(d, :w), ds) ?
    reduce(vcat, [hasproperty(d, :w) ? [(x=[xi], y=[yi], state_index=d.state_index, residual_scale=wi) for (xi, yi, wi) in zip(d.x, d.y, d.w)] :
                  [(x=d.x, y=d.y, state_index=d.state_index)] for d in ds]) : ds

function predict_candidate(name, label, model, datasets, u0)
    tmax = maximum(maximum(d.x) for d in datasets)
    sol = solve(ODEProblem(model, u0, (0.0, tmax), Float64[]), Tsit5(); reltol=1e-8, abstol=1e-8)
    sse = sum(sum(((sol(xi)[d.state_index] - yi) / wi)^2 for (xi, yi, wi) in zip(d.x, d.y, dweights(d))) for d in datasets)
    n = sum(length(d.x) for d in datasets)
    (name=name, label=label, model=model, pnames=String[], params=Float64[], bic=n*log(max(sse, 1e-12)/n), sse=sse, bounds=Tuple{Float64,Float64}[])
end

# Asymptotic standard errors from the least-squares Jacobian (a lower bound; residuals are autocorrelated).
function param_se(model, datasets, u0, p)
    k = length(p); k == 0 && return Float64[]
    tmax = maximum(maximum(d.x) for d in datasets)
    resid(q) = (sol = solve(ODEProblem(model, u0, (0.0, tmax), q), Tsit5(); reltol=1e-10, abstol=1e-10);
                reduce(vcat, [[(sol(xi)[d.state_index] - yi) / wi for (xi, yi, wi) in zip(d.x, d.y, dweights(d))] for d in datasets]))
    r0 = resid(p); n = length(r0); J = zeros(n, k)
    for j in 1:k
        h = 1e-5 * max(abs(p[j]), 1e-2); qp = copy(p); qp[j] += h; qm = copy(p); qm[j] -= h
        J[:, j] = (resid(qp) - resid(qm)) / (2h)
    end
    sqrt.(max.(diag(sum(abs2, r0) / max(n - k, 1) * pinv(J' * J)), 0.0))
end

fmt(v) = abs(v) >= 100 ? @sprintf("%.1f", v) : @sprintf("%.4g", v)
