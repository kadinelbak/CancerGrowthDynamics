module HybridRefitService

using CSV
using DataFrames
using Dates
using HTTP
using JSON3
using SHA

using ..IOUtils
using ..FitWorkflows
using ..StagedA2780Workflow

export canonical_refit_request,
       refit_request_id,
       validate_refit_request,
       run_conditional_refit!,
       serve_hybrid_refits

const REQUIRED_GROUPS = (
    "stage1_A2780Naive", "stage1_A2780cis",
    "stage2_A2780Naive", "stage2_A2780cis",
    "stage3", "stage4",
)
const STAGE4_COMPATIBLE = Dict(
    "A2780Naive" => "joint_ic_effect_hill_ramp_onset",
    "A2780cis" => "joint_ic_effect_two_population",
)
const JOBS = Dict{String,Dict{String,Any}}()
const JOB_LOCK = ReentrantLock()

_string(value) = String(value)
_now() = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"

function _selection(value)
    hasproperty(value, :model) || error("Each selection requires model")
    hasproperty(value, :pooling) || error("Each selection requires pooling")
    return (model = _string(value.model), pooling = _string(value.pooling))
end

"Return a stable, minimal representation used by both hashing and execution."
function canonical_refit_request(request)
    hasproperty(request, :selections) || error("Request requires selections")
    selections = request.selections
    pairs = Pair{Symbol,Any}[]
    for key in REQUIRED_GROUPS
        symbol = Symbol(key)
        hasproperty(selections, symbol) || error("Missing selection $(key)")
        push!(pairs, symbol => _selection(getproperty(selections, symbol)))
    end
    max_time = hasproperty(request, :max_time_per_fit) ? Float64(request.max_time_per_fit) : 12.0
    0.1 <= max_time <= 120.0 || error("max_time_per_fit must be between 0.1 and 120 seconds")
    return (schema_version = 1, selections = (; pairs...), max_time_per_fit = max_time)
end

refit_request_id(request) = bytes2hex(sha256(JSON3.write(canonical_refit_request(request))))[1:20]

function _canonical_artifacts(start)
    root = IOUtils.package_root(start)
    csv_root = joinpath(root, "outputs", "csv")
    return (
        root = root,
        stage1_ranking = joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_pooling_model_ranking.csv"),
        stage1_parameters = joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_pooling_parameter_estimates.csv"),
    )
end

function _selected_row(ranking, selection; cell_line = nothing)
    mask = (String.(ranking.model) .== selection.model) .&
           (String.(ranking.pooling_mode) .== selection.pooling)
    cell_line === nothing || (mask .&= String.(ranking.cell_line) .== cell_line)
    rows = ranking[mask, :]
    nrow(rows) == 1 || error("Expected one fitted row for $(selection.model) / $(selection.pooling)")
    row = first(rows)
    (:eligible_for_inheritance in propertynames(ranking) && !Bool(row.eligible_for_inheritance)) &&
        error("$(selection.model) / $(selection.pooling) is diagnostic and cannot be inherited")
    isfinite(Float64(row.bic)) || error("Selected fit has non-finite BIC")
    return row
end

function validate_refit_request(request; start = pwd())
    canonical = canonical_refit_request(request)
    artifacts = _canonical_artifacts(start)
    isfile(artifacts.stage1_ranking) || error("Stage 1 ranking artifact is missing")
    ranking = CSV.read(artifacts.stage1_ranking, DataFrame)
    _selected_row(ranking, canonical.selections.stage1_A2780Naive; cell_line = "A2780Naive")
    _selected_row(ranking, canonical.selections.stage1_A2780cis; cell_line = "A2780cis")
    stage2_compatible = all(
        canonical.selections[Symbol("stage2_$(lineage)")].model == model
        for (lineage, model) in STAGE4_COMPATIBLE
    )
    stage3_compatible = canonical.selections.stage3.model in (
        "lv_symmetric_competition", "lv_asymmetric_competition", "lv_asymmetric_competition_death",
    )
    return (request = canonical, stage2_compatible = stage2_compatible,
            stage3_compatible = stage3_compatible,
            stage4_compatible = stage2_compatible && stage3_compatible)
end

function _prepare_workspace(source_start, job_root)
    source_repo = IOUtils.find_repo_root(source_start)
    workspace = joinpath(job_root, "workspace")
    if !isdir(workspace)
        mkpath(workspace)
        cp(joinpath(source_repo, "Processed_Datasets"), joinpath(workspace, "Processed_Datasets"); force = true)
        mkpath(joinpath(workspace, "Modeling_Approaches", "02_mechanical_automatic_package", "outputs"))
    end
    return workspace
end

function _parameter_value(rows, name)
    selected = rows[String.(rows.parameter) .== name, :]
    return isempty(selected) ? NaN : Float64(first(selected.effective_value))
end

function _write_stage1_baseline!(workspace, source_start, selections)
    artifacts = _canonical_artifacts(source_start)
    ranking = CSV.read(artifacts.stage1_ranking, DataFrame)
    parameters = CSV.read(artifacts.stage1_parameters, DataFrame)
    baseline_rows = NamedTuple[]
    for cell_line in ("A2780Naive", "A2780cis")
        selection = selections[Symbol("stage1_$(cell_line)")]
        winner = _selected_row(ranking, selection; cell_line = cell_line)
        selected = parameters[
            (String.(parameters.cell_line) .== cell_line) .&
            (String.(parameters.model) .== selection.model) .&
            (String.(parameters.pooling_mode) .== selection.pooling), :]
        isempty(selected) && error("No Stage 1 parameters found for $(cell_line)")
        for density in sort(unique(String.(selected.density)))
            rows = selected[String.(selected.density) .== density, :]
            r = _parameter_value(rows, "r")
            K = _parameter_value(rows, "K")
            all(isfinite, (r, K)) || error("Stage 1 $(cell_line) is missing finite r/K at $(density)")
            shape_index = findfirst(name -> isfinite(_parameter_value(rows, name)), ["theta", "death_rate", "allee_threshold"])
            shape_name = shape_index === nothing ? "" : ["theta", "death_rate", "allee_threshold"][shape_index]
            push!(baseline_rows, (
                cell_line = cell_line, density = density, best_model = selection.model,
                pooling_mode = selection.pooling, r = r, K = K,
                shape_parameter = shape_name,
                shape_value = isempty(shape_name) ? NaN : _parameter_value(rows, shape_name),
                theta = _parameter_value(rows, "theta"), lag_time = _parameter_value(rows, "lag_time"),
                baranyi_q0 = _parameter_value(rows, "q0"), adaptation_rate = _parameter_value(rows, "adaptation_rate"),
                bic = Float64(winner.bic), ssr = Float64(winner.ssr), scaled_ssr = Float64(winner.scaled_ssr),
                inheritance_allowed = true, independent_bic_improvement = NaN,
            ))
        end
    end
    out = IOUtils.condition_output_dirs("monoculture_untreated"; start = workspace)
    CSV.write(joinpath(out.csv, "untreated_group_baselines.csv"), DataFrame(baseline_rows))
    return DataFrame(baseline_rows)
end

function _decode_and_fit(condition, workspace, max_time)
    if lowercase(get(ENV, "A2780_TOURNAMENT_SMOKE", "false")) == "true" && condition != "coculture_treated"
        output = IOUtils.condition_output_dirs(condition; start = workspace).csv
        ranking_name = condition == "monoculture_treated" ?
            "monoculture_treated_pooling_model_ranking.csv" :
            "coculture_untreated_pooling_model_ranking.csv"
        ranking_path = joinpath(output, ranking_name)
        isfile(ranking_path) && return (ranking = CSV.read(ranking_path, DataFrame), resumed = true)
    end
    decoded = StagedA2780Workflow.decode_a2780_condition(condition; start = workspace)
    out = IOUtils.condition_output_dirs(condition; start = workspace)
    CSV.write(joinpath(out.csv, "$(condition)_a2780_decoded.csv"), decoded)
    return FitWorkflows.run_condition_fit!(decoded, condition; start = workspace, max_time_per_fit = max_time)
end

function _write_stage2_selection!(workspace, selections, ranking)
    rows = NamedTuple[]
    for cell_line in ("A2780Naive", "A2780cis")
        selection = selections[Symbol("stage2_$(cell_line)")]
        winner = _selected_row(ranking, selection; cell_line = cell_line)
        push!(rows, (
            cell_line = cell_line, winning_model = selection.model,
            winning_pooling_mode = selection.pooling, winning_bic = Float64(winner.bic),
            best_independent_bic = NaN, independent_bic_improvement = NaN,
            inadequacy_delta = 10.0, inadequate_pooling = false, inheritance_allowed = true,
        ))
    end
    path = joinpath(IOUtils.condition_output_dirs("monoculture_treated"; start = workspace).csv,
                    "monoculture_treated_pooling_status.csv")
    CSV.write(path, DataFrame(rows))
end

function _values_from_params(text)
    matched = match(r"values = \[([^\]]+)\]", String(text))
    matched === nothing && error("Could not recover fitted parameter vector")
    return parse.(Float64, strip.(split(matched.captures[1], ",")))
end

function _write_stage3_selection!(workspace, selection, ranking)
    winner = _selected_row(ranking, selection)
    values = _values_from_params(winner.params)
    names = if selection.model == "lv_symmetric_competition"
        selection.pooling == "partial_5pct" ? ["alpha", "log_contrast_density"] : ["alpha"]
    elseif selection.model == "lv_asymmetric_competition"
        selection.pooling == "partial_5pct" ? ["alpha_sr", "alpha_rs", "log_contrast_density"] : ["alpha_sr", "alpha_rs"]
    elseif selection.model == "lv_asymmetric_competition_death"
        selection.pooling == "partial_5pct" ? ["alpha_sr", "alpha_rs", "death_sensitive", "death_resistant", "log_contrast_density"] : ["alpha_sr", "alpha_rs", "death_sensitive", "death_resistant"]
    else
        error("Stage 4 currently supports inherited LV competition models; $(selection.model) can be previewed but not inherited")
    end
    length(names) == length(values) || error("Selected Stage 3 parameter vector does not match its equation")
    baseline = DataFrame(
        model = [selection.model], pooling_mode = [selection.pooling],
        parameter_names = [join(names, ";")], parameter_values = [join(values, ";")],
        bic = [Float64(winner.bic)], ssr = [Float64(winner.ssr)], scaled_ssr = [Float64(winner.scaled_ssr)],
        monoculture_growth_inheritance = ["selected Stage 1 density-specific growth equations and parameters"],
    )
    path = joinpath(IOUtils.condition_output_dirs("coculture_untreated"; start = workspace).csv,
                    "coculture_untreated_joint_baseline.csv")
    CSV.write(path, baseline)
end

_json_value(value::AbstractFloat) = isfinite(value) ? value : nothing
_json_value(::Missing) = nothing
_json_value(value) = value

function _records(df; limit = nrow(df))
    return [Dict(String(name) => _json_value(row[Symbol(name)]) for name in names(df)) for row in eachrow(first(df, min(limit, nrow(df))))]
end

function _result_payload(id, canonical, workspace, stage2, stage3; stage4 = nothing, state = "completed", note = "")
    result = Dict{String,Any}(
        "job_id" => id, "status" => state, "completed_at" => _now(), "note" => note,
        "selection" => canonical, "workspace" => workspace,
        "stage2_ranking" => _records(sort(stage2.ranking, :bic); limit = 20),
        "stage3_ranking" => _records(sort(stage3.ranking, :bic); limit = 20),
        "provenance" => "conditionally refitted with MechanicalAutomaticModeling package functions",
    )
    stage4 === nothing || (result["stage4_ranking"] = _records(sort(stage4.ranking, :bic); limit = 20))
    return result
end

function _write_json(path, value)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, value)
    end
    return path
end

"Run the selected inheritance path in an isolated, content-addressed workspace."
function run_conditional_refit!(request; start = pwd(), cache_root = nothing, force = false)
    validation = validate_refit_request(request; start = start)
    canonical = validation.request
    id = refit_request_id(canonical)
    root = something(cache_root, joinpath(IOUtils.package_root(start), "outputs", "refit_cache"))
    job_root = joinpath(root, id)
    result_path = joinpath(job_root, "result.json")
    if isfile(result_path) && !force
        return JSON3.read(read(result_path, String))
    end
    workspace = _prepare_workspace(start, job_root)
    _write_json(joinpath(job_root, "request.json"), canonical)
    _write_stage1_baseline!(workspace, start, canonical.selections)
    stage2 = _decode_and_fit("monoculture_treated", workspace, canonical.max_time_per_fit)
    _write_stage2_selection!(workspace, canonical.selections, stage2.ranking)
    stage3 = _decode_and_fit("coculture_untreated", workspace, canonical.max_time_per_fit)
    validation.stage3_compatible && _write_stage3_selection!(workspace, canonical.selections.stage3, stage3.ranking)
    if !validation.stage4_compatible
        reasons = String[]
        validation.stage2_compatible || push!(reasons, "Stage 4 requires delayed Hill-ramp A2780Naive and sensitive/tolerant A2780cis Stage 2 states")
        validation.stage3_compatible || push!(reasons, "Stage 4 currently inherits one of the three LV competition equations")
        note = "Stages 1-3 were conditionally refitted. " * join(reasons, "; ") * "."
        result = _result_payload(id, canonical, workspace, stage2, stage3; state = "partial", note = note)
        _write_json(result_path, result)
        return result
    end
    # Make the linked seed use the exact selected pooling rather than another eligible pooling.
    ranking_path = joinpath(IOUtils.condition_output_dirs("monoculture_treated"; start = workspace).csv,
                            "monoculture_treated_pooling_model_ranking.csv")
    stage2_ranking = CSV.read(ranking_path, DataFrame)
    for cell_line in ("A2780Naive", "A2780cis")
        selected = canonical.selections[Symbol("stage2_$(cell_line)")]
        mask = (String.(stage2_ranking.cell_line) .== cell_line) .&
               (String.(stage2_ranking.model) .== selected.model)
        stage2_ranking.eligible_for_inheritance[mask] .= String.(stage2_ranking.pooling_mode[mask]) .== selected.pooling
    end
    CSV.write(ranking_path, stage2_ranking)
    stage4 = _decode_and_fit("coculture_treated", workspace, canonical.max_time_per_fit)
    _selected_row(stage4.ranking, canonical.selections.stage4)
    result = _result_payload(id, canonical, workspace, stage2, stage3; stage4 = stage4)
    _write_json(result_path, result)
    return result
end

function _set_job!(id; kwargs...)
    lock(JOB_LOCK) do
        job = get!(JOBS, id, Dict{String,Any}("job_id" => id))
        for (key, value) in kwargs
            job[String(key)] = value
        end
    end
end

function _job(id, cache_root)
    lock(JOB_LOCK) do
        haskey(JOBS, id) && return copy(JOBS[id])
    end
    result_path = joinpath(cache_root, id, "result.json")
    if isfile(result_path)
        result = JSON3.read(read(result_path, String))
        return Dict{String,Any}("job_id" => id, "status" => String(result.status), "cached" => true, "result" => result)
    end
    return nothing
end

_cors() = ["Access-Control-Allow-Origin" => "*", "Access-Control-Allow-Headers" => "Content-Type", "Access-Control-Allow-Methods" => "GET, POST, OPTIONS"]
_json_response(status, value) = HTTP.Response(status, _cors(), JSON3.write(value))

function serve_hybrid_refits(; host = "127.0.0.1", port = 8766, start = pwd())
    cache_root = joinpath(IOUtils.package_root(start), "outputs", "refit_cache")
    mkpath(cache_root)
    handler = function(request)
        request.method == "OPTIONS" && return HTTP.Response(204, _cors())
        path = HTTP.URI(request.target).path
        if request.method == "GET" && path == "/health"
            return _json_response(200, (status = "ok", service = "A2780 conditional refit", api_version = 1))
        elseif request.method == "POST" && path == "/api/refits"
            try
                body = JSON3.read(String(request.body))
                canonical = canonical_refit_request(body)
                validation = validate_refit_request(canonical; start = start)
                id = refit_request_id(canonical)
                cached = isfile(joinpath(cache_root, id, "result.json"))
                if !cached && _job(id, cache_root) === nothing
                    _set_job!(id; status = "queued", submitted_at = _now(), stage4_compatible = validation.stage4_compatible)
                    Threads.@spawn begin
                        _set_job!(id; status = "running", started_at = _now())
                        try
                            result = run_conditional_refit!(canonical; start = start, cache_root = cache_root)
                            _set_job!(id; status = String(result["status"]), completed_at = _now(), result = result)
                        catch error_value
                            _set_job!(id; status = "failed", completed_at = _now(), error = sprint(showerror, error_value))
                        end
                    end
                end
                job = _job(id, cache_root)
                return _json_response(cached ? 200 : 202, merge(Dict("cached" => cached), job))
            catch error_value
                return _json_response(400, (status = "invalid", error = sprint(showerror, error_value)))
            end
        elseif request.method == "GET" && startswith(path, "/api/refits/")
            id = split(path, '/')[end]
            job = _job(id, cache_root)
            job === nothing && return _json_response(404, (status = "missing", job_id = id))
            return _json_response(200, job)
        end
        return _json_response(404, (status = "not_found", path = path))
    end
    println("A2780 conditional-refit service: http://$(host):$(port)")
    return HTTP.serve(handler, host, port)
end

end
