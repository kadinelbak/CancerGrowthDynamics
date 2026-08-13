module IOUtils

using CSV
using DataFrames
using Dates

export find_repo_root,
       package_root,
       condition_processed_dir,
       condition_output_dirs,
       decode_condition_dataframe,
       write_manifest_row,
       sanitize_name

const CONDITION_DIRS = Dict(
    "monoculture_untreated" => joinpath("Processed_Datasets", "Untreated MonoCulture"),
    "monoculture_treated" => joinpath("Processed_Datasets", "Treated MonoCulture"),
    "coculture_untreated" => joinpath("Processed_Datasets", "Untreated CoCulture"),
    "coculture_treated" => joinpath("Processed_Datasets", "Treated Coculture"),
)

# The historical splitter wrote the 1.47 uM source into folders named IC25.
# Apply the later user-confirmed biological labels while leaving source files intact:
# legacy IC25 -> 1.47 uM (IC75), legacy IC75 -> 0.67 uM (IC25).
const TREATED_MONOCULTURE_IC50_UM = 1.0
const TREATED_MONOCULTURE_IC_DOSE_MAP = Dict(
    25 => 1.47,
    50 => 1.0,
    75 => 0.67,
)

function find_repo_root(start::AbstractString = pwd())
    current = abspath(start)
    for _ in 1:10
        if isdir(joinpath(current, "Processed_Datasets")) && isdir(joinpath(current, "Modeling_Approaches"))
            return current
        end
        if isdir(joinpath(current, "Modeling_Approaches")) && isdir(joinpath(current, "Processing Files"))
            return current
        end
        parent = dirname(current)
        parent == current && break
        current = parent
    end
    error("Repository root not found from: $(start)")
end

function package_root(start::AbstractString = pwd())
    repo = find_repo_root(start)
    return joinpath(repo, "Modeling_Approaches", "02_mechanical_automatic_package")
end

sanitize_name(s::AbstractString) = replace(lowercase(strip(s)), r"[^a-z0-9]+" => "_")

function condition_processed_dir(condition::AbstractString; start::AbstractString = pwd())
    haskey(CONDITION_DIRS, condition) || error("Unknown condition: $(condition)")
    return joinpath(find_repo_root(start), CONDITION_DIRS[condition])
end

function condition_output_dirs(condition::AbstractString; start::AbstractString = pwd())
    root = package_root(start)
    base = joinpath(root, "outputs")
    csv_dir = joinpath(base, "csv", condition)
    img_dir = joinpath(base, "images", condition)
    metric_dir = joinpath(base, "metrics", condition)
    mkpath(csv_dir)
    mkpath(img_dir)
    mkpath(metric_dir)
    return (csv = csv_dir, images = img_dir, metrics = metric_dir)
end

function _coerce_float(v)
    (v === missing || v === nothing) && return missing
    txt = strip(string(v))
    isempty(txt) && return missing
    return something(tryparse(Float64, txt), missing)
end

function _normalize_columns(df::DataFrame)
    rename!(df, Dict(c => Symbol(sanitize_name(String(c))) for c in names(df)))
    return df
end

function _first_existing_column(cols, candidates::Vector{Symbol})
    colset = Set(Symbol.(cols))
    for c in candidates
        c in colset && return c
    end
    return nothing
end

function _replicate_from_values(vals)
    labels = [strip(string(v)) for v in vals]
    idx = Dict{String,Int}()
    next_id = 1
    reps = Vector{Int}(undef, length(labels))
    for i in eachindex(labels)
        key = labels[i]
        if isempty(key)
            reps[i] = 1
            continue
        end
        if !haskey(idx, key)
            idx[key] = next_id
            next_id += 1
        end
        reps[i] = idx[key]
    end
    return reps
end

function _infer_metadata_from_path(f::AbstractString)
    lf = lowercase(f)
    parts = splitpath(f)

    density = ""
    for p in parts
        pl = lowercase(p)
        if occursin(r"^\d+k$", pl)
            density = p
            break
        end
    end

    cell_line = ""
    if occursin("a2780cis", lf)
        cell_line = "A2780cis"
    elseif occursin("a2780naive", lf)
        cell_line = "A2780Naive"
    elseif occursin("a2780", lf) && occursin("coculture", lf)
        cell_line = "A2780Naive"
    elseif occursin("a2780", lf)
        cell_line = "A2780"
    end

    dose = 0.0
    m_ic = match(r"ic(\d+)"i, lf)
    if m_ic !== nothing
        ic_level = parse(Int, m_ic.captures[1])
        if occursin("treated monoculture", lf)
            dose = get(TREATED_MONOCULTURE_IC_DOSE_MAP, ic_level, TREATED_MONOCULTURE_IC50_UM)
        elseif occursin("ic50", lf) || occursin("1um", lf)
            dose = 1.0
        else
            dose = Float64(ic_level)
        end
    else
        m_um = match(r"([0-9]+(?:\.[0-9]+)?)um"i, lf)
        if m_um !== nothing
            dose = parse(Float64, m_um.captures[1])
        end
    end

    mix = ""
    m_mix = match(r"(\d{2})-(\d{2})", lf)
    if m_mix !== nothing
        mix = "$(m_mix.captures[1])-$(m_mix.captures[2])"
    end

    return (density = density, cell_line = cell_line, dose = dose, mix = mix)
end

function decode_condition_dataframe(condition::AbstractString; start::AbstractString = pwd())
    pdir = condition_processed_dir(condition; start)
    files = String[]
    for (root, _, names_in_dir) in walkdir(pdir)
        for name in names_in_dir
            endswith(lowercase(name), ".csv") || continue
            push!(files, joinpath(root, name))
        end
    end
    sort!(files)
    isempty(files) && error("No CSV files found in $(pdir)")

    parts = DataFrame[]
    for f in files
        lf = lowercase(f)
        if endswith(lf, "_day_averages.csv")
            paired_sample = replace(f, r"(?i)_day_averages\.csv$" => "_sample_averages.csv")
            if isfile(paired_sample)
                continue
            end
        end

        df = CSV.read(f, DataFrame)
        _normalize_columns(df)
        inferred = _infer_metadata_from_path(f)

        cols = names(df)
        colset = Set(Symbol.(cols))
        tcol = _first_existing_column(cols, [:day, :time, :t, :days])
        ycol = _first_existing_column(cols, [:mean_cells, :cell_count, :count, :cells, :total_cells, :mean_count, :day_mean_value, :mean_value])

        if tcol === nothing || ycol === nothing
            continue
        end

        out = DataFrame(
            time = [_coerce_float(v) for v in df[!, tcol]],
            count = [_coerce_float(v) for v in df[!, ycol]],
            source_file = fill(basename(f), nrow(df)),
            condition = fill(condition, nrow(df)),
        )

        if :density in colset
            out.density = df[!, :density]
        elseif :seed_density in colset
            out.density = df[!, :seed_density]
        else
            out.density = fill(inferred.density, nrow(df))
        end

        if :cell_line in colset
            out.cell_line = df[!, :cell_line]
        elseif :cellline in colset
            out.cell_line = df[!, :cellline]
        else
            out.cell_line = fill(inferred.cell_line, nrow(df))
        end

        if :dose in colset
            out.dose = [_coerce_float(v) for v in df[!, :dose]]
        elseif :dose_um in colset
            out.dose = [_coerce_float(v) for v in df[!, :dose_um]]
        else
            out.dose = fill(inferred.dose, nrow(df))
        end

        is_coculture = occursin("coculture", lowercase(condition))
        if :mix in colset
            out.mix = df[!, :mix]
        elseif :mix_label in colset
            out.mix = df[!, :mix_label]
        elseif :well in colset && is_coculture
            out.mix = fill(inferred.mix, nrow(df))
        else
            out.mix = fill(inferred.mix, nrow(df))
        end

        if :replicate in colset
            reps = [_coerce_float(v) for v in df[!, :replicate]]
            out.replicate = [ismissing(v) ? 1 : max(1, Int(round(v))) for v in reps]
        elseif :well in colset
            out.replicate = _replicate_from_values(df[!, :well])
        elseif :mix in colset
            out.replicate = _replicate_from_values(df[!, :mix])
        elseif :mix_label in colset
            out.replicate = _replicate_from_values(df[!, :mix_label])
        else
            out.replicate = fill(1, nrow(df))
        end

        filter!(row -> !ismissing(row.time) && !ismissing(row.count), out)
        isempty(out) || push!(parts, out)
    end

    isempty(parts) && error("Could not decode any usable rows for condition: $(condition)")
    decoded = vcat(parts...)
    sort!(decoded, [:time, :source_file])
    return decoded
end

function write_manifest_row(; condition::AbstractString, step::AbstractString, outputs::Vector{String}, start::AbstractString = pwd())
    root = package_root(start)
    manifest_dir = joinpath(root, "outputs", "manifests")
    mkpath(manifest_dir)
    manifest_file = joinpath(manifest_dir, "run_manifest.csv")

    row = DataFrame(
        timestamp_utc = [Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")],
        condition = [condition],
        step = [step],
        output_count = [length(outputs)],
        outputs = [join(outputs, ";")],
    )

    if isfile(manifest_file)
        CSV.write(manifest_file, row; append = true)
    else
        CSV.write(manifest_file, row)
    end
    return manifest_file
end

end
