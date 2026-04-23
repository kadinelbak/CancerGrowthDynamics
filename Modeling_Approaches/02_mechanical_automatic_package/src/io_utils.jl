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

function find_repo_root(start::AbstractString = pwd())
    current = abspath(start)
    for _ in 1:10
        if isdir(joinpath(current, "Processed_Datasets")) && isdir(joinpath(current, "Modeling_Approaches"))
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
    if v === missing || v === nothing
        return missing
    end
    txt = strip(string(v))
    isempty(txt) && return missing
    x = tryparse(Float64, txt)
    x === nothing ? missing : x
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
        df = CSV.read(f, DataFrame)
        _normalize_columns(df)

        cols = names(df)
        tcol = _first_existing_column(cols, [:day, :time, :t, :days])
        ycol = _first_existing_column(cols, [:mean_cells, :cell_count, :count, :cells, :total_cells, :mean_count])

        if tcol === nothing || ycol === nothing
            continue
        end

        out = DataFrame(
            time = [_coerce_float(v) for v in df[!, tcol]],
            count = [_coerce_float(v) for v in df[!, ycol]],
            source_file = fill(basename(f), nrow(df)),
            condition = fill(condition, nrow(df)),
        )

        if :density in cols
            out.density = df[!, :density]
        elseif :seed_density in cols
            out.density = df[!, :seed_density]
        else
            out.density = fill("", nrow(df))
        end

        if :cell_line in cols
            out.cell_line = df[!, :cell_line]
        elseif :cellline in cols
            out.cell_line = df[!, :cellline]
        else
            out.cell_line = fill("", nrow(df))
        end

        if :dose in cols
            out.dose = [_coerce_float(v) for v in df[!, :dose]]
        elseif :dose_um in cols
            out.dose = [_coerce_float(v) for v in df[!, :dose_um]]
        else
            out.dose = fill(0.0, nrow(df))
        end

        if :mix in cols
            out.mix = df[!, :mix]
        elseif :mix_label in cols
            out.mix = df[!, :mix_label]
        else
            out.mix = fill("", nrow(df))
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
