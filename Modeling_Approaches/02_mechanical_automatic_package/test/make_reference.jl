# Saves the results of the current run as the regression reference (test/reference_results.toml):
# every candidate's BIC and the carried model's parameters per group, the first bootstrap
# replicates' BICs, and the control outcome tables. Run once after a full run you accept:
#   julia --project=. test/make_reference.jl
using Serialization, TOML, CSV, DataFrames
const PK = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PK, "zoo_common.jl")); include(joinpath(PK, "zoo_models.jl"))
state = deserialize(joinpath(PK, "outputs", "zoo_fits", "zoo_state.jls"))
# configuration fingerprint: a reference only applies to the configuration and data it was made with
include(joinpath(@__DIR__, "regression_fingerprint.jl"))
ref = Dict{String,Any}("made" => string(now()), "groups" => Dict{String,Any}(), "fingerprint" => config_fingerprint())
for g in state.groups
    r = state.results[g.id].res; nb = min(3, size(r.bic, 1))
    ref["groups"][g.id] = Dict("carried" => g.cands[1].name,
                               "bic" => Dict(c.name => c.bic for c in g.cands),
                               "params" => Dict(n => v for (n, v) in zip(g.cands[1].pnames, g.cands[1].params)),
                               "boot_bic" => [Dict(c.name => r.bic[b, j] for (j, c) in enumerate(r.cands)) for b in 1:nb])
end
const LONG = round(Int, TOML.parsefile(joinpath(PK, "control_config.toml"))["schedule"]["long_days"])
for (name, f) in (("control14", "outputs/control/control_outcomes.csv"), ("control_long", "outputs/control_$(LONG)d/fasting_outcomes.csv"))
    isfile(joinpath(PK, f)) || continue
    df = CSV.read(joinpath(PK, f), DataFrame)
    ref[name] = [Dict(string(k) => (v isa Number ? Float64(v) : string(v)) for (k, v) in pairs(row) if !ismissing(v)) for row in eachrow(df)]
end
open(io -> TOML.print(io, ref), joinpath(@__DIR__, "reference_results.toml"), "w")
println("Wrote test/reference_results.toml ($(length(state.groups)) groups)")
