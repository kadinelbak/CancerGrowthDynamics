using Pkg
Pkg.activate(pwd())
using CSV, DataFrames, Statistics
include(joinpath(pwd(), "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

condition = "monoculture_untreated"
out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs(condition; start=joinpath(pwd(), "notebooks"))

decoded = MechanicalAutomaticModeling.IOUtils.decode_condition_dataframe(condition; start=joinpath(pwd(), "notebooks"))
sort!(decoded, [:cell_line, :density, :dose, :replicate, :time])
CSV.write(joinpath(out.csv, "$(condition)_automatic_decoded.csv"), decoded)

group_cols = Symbol[]
:cell_line in names(decoded) && push!(group_cols, :cell_line)
:density in names(decoded) && push!(group_cols, :density)
:dose in names(decoded) && push!(group_cols, :dose)
coverage = combine(groupby(decoded, group_cols), nrow => :n_rows, :replicate => (x -> length(unique(x))) => :n_replicates)
CSV.write(joinpath(out.metrics, "$(condition)_decoded_group_coverage.csv"), coverage)

fit_artifacts = MechanicalAutomaticModeling.FitWorkflows.run_condition_fit!(decoded, condition; start=joinpath(pwd(), "notebooks"))

sanitize_tag(s) = replace(lowercase(string(s)), r"[^a-z0-9]+" => "_")
parse_params(params_txt) = [parse(Float64, m.match) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", String(params_txt))]
function model_prediction(model_name::String, params::Vector{Float64}, t::Vector{Float64}, n0::Float64)
    t0 = first(t)
    if model_name == "logistic_growth" && length(params) >= 2
        r, K = params[1], max(params[2], 1e-8)
        return [K / (1 + ((K - n0) / max(n0, 1e-8)) * exp(-r * (ti - t0))) for ti in t]
    elseif model_name == "gompertz_growth" && length(params) >= 2
        r, K = params[1], max(params[2], 1e-8)
        return [K * exp(log(max(n0, 1e-8) / K) * exp(-r * (ti - t0))) for ti in t]
    else
        return fill(NaN, length(t))
    end
end

overlay_img_dir = joinpath(out.images, "figures")
overlay_csv_dir = joinpath(out.csv, "figures")
mkpath(overlay_img_dir)
mkpath(overlay_csv_dir)

saved_paths = String[]
frames = DataFrame[]
for grp in groupby(decoded, group_cols)
    key_map = Dict{Symbol,Any}(c => first(grp[!, c]) for c in group_cols)
    group_rank = fit_artifacts.ranking
    for c in group_cols
        if c in names(group_rank)
            group_rank = filter(r -> r[c] == key_map[c], group_rank)
        end
    end
    nrow(group_rank) == 0 && continue

    best = first(sort(group_rank, :bic), 1)
    model = String(best.model[1])
    pvec = parse_params(best.params[1])
    obs_mean = combine(groupby(grp, :time), :count => mean => :observed)
    sort!(obs_mean, :time)
    pred = model_prediction(model, pvec, Float64.(obs_mean.time), max(first(obs_mean.observed), 1e-8))

    overlay_df = DataFrame(time=obs_mean.time, observed=obs_mean.observed, prediction=pred, best_model=fill(model, nrow(obs_mean)))
    for c in group_cols
        overlay_df[!, c] = fill(string(key_map[c]), nrow(obs_mean))
    end

    tag = join(["$(c)_$(sanitize_tag(key_map[c]))" for c in group_cols], "__")
    csv_path = joinpath(overlay_csv_dir, "$(condition)_$(tag)_best_overlay.csv")
    CSV.write(csv_path, overlay_df)
    push!(saved_paths, csv_path)
    push!(frames, overlay_df)
end

if !isempty(frames)
    all_path = joinpath(overlay_csv_dir, "$(condition)_all_groups_best_overlays.csv")
    CSV.write(all_path, vcat(frames...))
    push!(saved_paths, all_path)
end

println("decoded_rows=$(nrow(decoded))")
println("n_groups=$(nrow(coverage))")
println("fit_rows=$(nrow(fit_artifacts.ranking))")
println("overlay_exports=$(length(saved_paths))")
println("sample_overlay_file=" * (isempty(saved_paths) ? "none" : basename(saved_paths[1])))
