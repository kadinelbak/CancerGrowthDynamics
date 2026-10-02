# Stage 0: build every 20k and 30k A2780 / A2780cis trajectory in the repo, from the rawest source
# available, and apply documented quality control. No GrowthParameterEstimation.jl.
#
# Sources (audit 2026-09-28; see outputs/data/report.html)
#   raw per-image files (Processed_Datasets): untreated + treated monoculture, intermittent,
#     untreated (ce0) + treated (ce1) co-culture -> per-well day means -> day mean/sd/sem.
#     Co-culture wells are assigned by plate column (Well Layouts: cols 1-3 = 20k, 4-6 = 30k;
#     75/25, 50/50, 25/75 by column). The legacy co-culture Averages split by tile and mix densities.
#   LR workbooks (no raw images): low-resource runs 1 and 2, 30k sheets only.
#   New Datasets monoculture CSVs and ce0/ce1 workbooks duplicate the raw data and are not loaded.
#   Both 20k and 30k seeding are used (user, 2026-09-29); seeding density is a covariate.
#   Not used: 40k wells, TykNu, New Datasets/Validation (LR + cisplatin), day 0 (a seeding value).
#
# Splits
#   :train                   used for fitting
#   :validation_treated_co   ce1 treated co-culture (held out)
# Roles inside :train
#   :drug         monoculture (+ intermittent) data: backbone + drug response
#   :competition  untreated co-culture: only the competition coefficients

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PD_DIR = joinpath(REPO_ROOT, "Processed_Datasets")
const LR_CSV = joinpath(REPO_ROOT, "Modeling_Approaches", "02_mechanical_automatic_package", "outputs",
                        "low_resource_workbook_trajectories.csv")
const INTERMITTENT_DIR = joinpath(PD_DIR, "Intermittent Data")

"Piecewise-constant dose schedule: dose `doses[i]` applies from `starts[i]` onward."
struct Schedule{T <: Real}
    starts::Vector{Float64}
    doses::Vector{T}
end
Schedule(d::Real) = Schedule([0.0], [Float64(d)])
dose_at(s::Schedule, t) = s.doses[searchsortedlast(s.starts, t) |> i -> max(i, 1)]
describe(s::Schedule) = length(s.doses) == 1 ? "const $(s.doses[1]) uM" :
    join(["d$(Int(a))+: $(b) uM" for (a, b) in zip(s.starts, s.doses)], ", ")

Base.@kwdef struct Trajectory
    id::String
    source::Symbol          # :mono, :lr, :ce0, :ce1, :intermittent
    plate::String           # experiment/plate, for plate-wide artefact detection
    split::Symbol
    role::Symbol
    line::Symbol            # :S (A2780) or :R (A2780cis)
    resource::Symbol        # :normal or :low
    culture::Symbol         # :mono or :co
    group::String           # well id; co-culture partners share it
    s_frac::Float64         # sensitive fraction of the seeded well
    seeding::Int            # 20 or 30 (k cells)
    units::Symbol           # :count or :area
    schedule::Schedule
    days::Vector{Float64}
    mean::Vector{Float64}
    sd::Vector{Float64}
    sem::Vector{Float64}
end

treated(tr::Trajectory) = any(>(0), tr.schedule.doses)

function _clean_names!(df)
    rename!(df, [strip(replace(String(n), '﻿' => "")) for n in names(df)])
    return df
end

"Every excluded point with the rule that removed it (for plotting/reporting)."
const EXCLUDED_POINTS = DataFrame(traj_id = String[], day = Float64[], mean = Float64[], rule = String[])
"Every excluded well-day (raw sources only)."
const EXCLUDED_WELLS = DataFrame(traj_id = String[], day = Float64[], well = String[], value = Float64[], rule = String[])

# ---------------- raw per-image files -> per-well day means ----------------
const IMG_RE = r"Day(\d+).*Tile-(\d+)_([A-D])(\d)\.tif"i

"Raw image file -> Dict(day => Dict(well => mean over tiles)), optionally restricted to plate columns."
function raw_well_days(path; columns = nothing)
    df = CSV.read(path, DataFrame; header = 1, silencewarnings = true)
    out = Dict{Float64, Dict{String, Vector{Float64}}}()
    for row in eachrow(df)
        img, val = row[1], row[2]
        (ismissing(img) || ismissing(val)) && continue
        m = match(IMG_RE, String(img))
        m === nothing && continue
        col = parse(Int, m[4])
        columns === nothing || col in columns || continue
        v = Float64(val)
        isfinite(v) || continue
        push!(get!(get!(out, parse(Float64, m[1]), Dict{String, Vector{Float64}}()), m[3] * m[4], Float64[]), v)
    end
    return Dict(d => Dict(w => mean(v) for (w, v) in ws) for (d, ws) in out)
end

const MONO_DOSE = Dict("Untreated" => 0.0, "IC25" => 1.47, "IC50" => 1.0, "IC75" => 0.62)  # folder -> uM (matches New Datasets file names)

"""
Raw sources: (id, source, plate, line, culture, s_frac, schedule, role, split, path, columns, group, seeding).
Monoculture and intermittent files are already split by density (well rows A/C = 20k, B/D = 30k);
co-culture flat files are split here by plate column (1-3 = 20k, 4-6 = 30k). 20k and 30k wells of the
same experiment share a plate, so plate-wide artefact detection pools them.
"""
function raw_specs()
    specs = NamedTuple[]
    for seeding in (20, 30), (cell, line) in (("A2780Naive", :S), ("A2780cis", :R))
        nm = cell == "A2780Naive" ? "A2780" : "A2780cis"
        k = "$(seeding)k"
        sf = line == :S ? 1.0 : 0.0
        push!(specs, (id = "mono_$(nm)_0uM_$k", source = :mono, plate = "mono untreated", line, culture = :mono,
            s_frac = sf, schedule = Schedule(0.0), role = :drug, split = :train, seeding,
            path = joinpath(PD_DIR, "Untreated MonoCulture", k, "$(cell).csv"), columns = nothing, group = "mono_$(nm)_0uM_$k"))
        for ic in ("IC25", "IC50", "IC75")
            d = MONO_DOSE[ic]
            plate = ic == "IC25" ? "mono 1.47uM ($(nm))" : "mono 1 + 0.62uM"
            push!(specs, (id = "mono_$(nm)_$(d)uM_$k", source = :mono, plate, line, culture = :mono,
                s_frac = sf, schedule = Schedule(d), role = :drug, split = :train, seeding,
                path = joinpath(PD_DIR, "Treated MonoCulture", k, ic, "$(cell).csv"), columns = nothing, group = "mono_$(nm)_$(d)uM_$k"))
        end
        for tag in ("T", "UT")
            # Understanding.md + Well Layouts: T = 1 uM days 1-7 then off; UT = off then 1 uM days 8-14.
            sched = tag == "T" ? Schedule([0.0, 7.0], [1.0, 0.0]) : Schedule([0.0, 7.0], [0.0, 1.0])
            id = "int_$(nm)_$(tag == "T" ? "on-off" : "off-on")_$k"
            push!(specs, (id, source = :intermittent, plate = "intermittent", line, culture = :mono,
                s_frac = sf, schedule = sched, role = :drug, split = :train, seeding,
                path = joinpath(INTERMITTENT_DIR, "$(k)_seeding_density", "$(nm)$(tag).csv"), columns = nothing, group = id))
        end
    end
    cols = Dict(20 => Dict("75-25" => 1, "50-50" => 2, "25-75" => 3),        # Well Layouts
                30 => Dict("75-25" => 4, "50-50" => 5, "25-75" => 6))
    for seeding in (20, 30), (ce, folder, pattern) in (("ce0", "Untreated CoCulture", r"^measure_(\d+)-(\d+)(A2780cis|A2780)CoUntreat"),
                                                       ("ce1", "Treated Coculture", r"^measure_(\d+)-(\d+)(A2780cis|A2780)IC50tre"))
        for f in sort(readdir(joinpath(PD_DIR, folder)))
            m = match(pattern, f)
            m === nothing && continue
            a, b, cell = parse(Int, m[1]), parse(Int, m[2]), m[3]
            line = cell == "A2780cis" ? :R : :S
            treated_ = ce == "ce1"
            group = "$(ce)_$(a)S$(b)R_$(seeding)k"
            push!(specs, (id = "$(group)_$(line)", source = Symbol(ce), plate = ce == "ce0" ? "co-culture untreated" : "co-culture treated",
                line, culture = :co, s_frac = a / (a + b), schedule = Schedule(treated_ ? 1.0 : 0.0),
                role = treated_ ? :validation : :competition, split = treated_ ? :validation_treated_co : :train, seeding,
                path = joinpath(PD_DIR, folder, f), columns = (cols[seeding]["$(a)-$(b)"],), group))
        end
    end
    return specs
end

"""
Well-level QC (temporal): a well-day is dropped when that well jumps > 2-fold against the mean of
its own previous and next day (a spike or a dip), while the other wells on that day move < 1.5-fold
the same way. Wells that are consistently high or low (seeding differences) are kept.
"""
function qc_wells(id, wd)
    days = sort(collect(keys(wd)))
    clean = Dict(d => copy(wd[d]) for d in days)
    for (k, d) in enumerate(days)
        (k == 1 || k == length(days)) && continue
        dp, dn = days[k - 1], days[k + 1]
        dip(w) = (haskey(wd[dp], w) && haskey(wd[dn], w) && min(wd[d][w], wd[dp][w], wd[dn][w]) > 0) ?
                 log(wd[d][w]) - (log(wd[dp][w]) + log(wd[dn][w])) / 2 : NaN
        dips = Dict(w => dip(w) for w in keys(wd[d]))
        for (w, x) in dips
            isnan(x) && continue
            others = [y for (v, y) in dips if v != w && !isnan(y)]
            isempty(others) && continue
            if abs(x) > log(2.0) && abs(mean(others)) < log(1.5)
                delete!(clean[d], w)
                push!(EXCLUDED_WELLS, (id, d, w, wd[d][w], x > 0 ? "well spike (>2x vs own neighbours)" : "well dip (>2x vs own neighbours)"))
            end
        end
    end
    return clean
end

function trajectory_from_wells(spec, wd)
    days = sort([d for d in keys(wd) if d >= 1 && !isempty(wd[d])])
    mu = [mean(values(wd[d])) for d in days]
    sd = [length(wd[d]) > 1 ? std(collect(values(wd[d]))) : 0.0 for d in days]
    sem = [length(wd[d]) > 1 ? sd[i] / sqrt(length(wd[d])) : 0.0 for (i, d) in enumerate(days)]
    return Trajectory(; id = spec.id, source = spec.source, plate = spec.plate, split = spec.split, role = spec.role,
        line = spec.line, resource = :normal, culture = spec.culture, group = spec.group, s_frac = spec.s_frac,
        seeding = spec.seeding, units = spec.source == :intermittent ? :area : :count, schedule = spec.schedule,
        days, mean = mu, sd, sem)
end

function load_raw()
    out = Trajectory[]
    for spec in raw_specs()
        wd = raw_well_days(spec.path; columns = spec.columns)
        push!(out, trajectory_from_wells(spec, qc_wells(spec.id, wd)))
    end
    return out
end

# ---------------- low-resource workbooks (no raw images available) ----------------
function load_lr()
    raw = CSV.read(LR_CSV, DataFrame)
    out = Trajectory[]
    for g in groupby(raw, [:workbook, :sheet])
        wb, sheet = String(g.workbook[1]), String(g.sheet[1])
        startswith(wb, "ALL A2780") || continue
        m = match(r"^(\d) - (mono )?(A2780cis|A2780) (\d+)k(?: (\d)-(\d))?$", sheet)
        m === nothing && error("unparsed LR sheet: $sheet")
        run, mono, cell, seeding = m[1], m[2] !== nothing, m[3], parse(Int, m[4])
        seeding in (20, 30) || continue
        series(name) = (s = sort(g[g.series .== name, [:day, :value]], :day); Dict(zip(Float64.(s.day), Float64.(s.value))))
        mu_d, sd_d, sem_d = series("Mean Cells"), series("SD Cells"), series("SEM Cells")
        days = sort([d for d in keys(mu_d) if d >= 1])          # day 0 is the seeding value, not a measurement
        line = cell == "A2780cis" ? :R : :S
        if mono
            s_frac, group = line == :S ? 1.0 : 0.0, "lr$(run)_mono_$(cell)_$(seeding)k"
        else
            a, b = parse(Int, m[5]), parse(Int, m[6])
            s_frac, group = a / (a + b), "lr$(run)_$(a)-$(b)_$(seeding)k"
        end
        push!(out, Trajectory(; id = "lr$(run)_$(cell)_$(mono ? "mono" : "$(m[5])-$(m[6])")_$(seeding)k", source = :lr,
            plate = "LR run $(run)", split = :train, role = mono ? :drug : :competition, line, resource = :low,
            culture = mono ? :mono : :co, group, s_frac, seeding, units = :count, schedule = Schedule(0.0),
            days, mean = [mu_d[d] for d in days], sd = [get(sd_d, d, 0.0) for d in days],
            sem = [get(sem_d, d, 0.0) for d in days]))
    end
    return out
end

# ---------------- point-level QC ----------------
function keep_points(tr::Trajectory, keep::AbstractVector{Bool}, rule)
    for i in findall(!, keep)
        push!(EXCLUDED_POINTS, (tr.id, tr.days[i], tr.mean[i], rule))
    end
    return Trajectory(; id = tr.id, source = tr.source, plate = tr.plate, split = tr.split, role = tr.role, line = tr.line,
        resource = tr.resource, culture = tr.culture, group = tr.group, s_frac = tr.s_frac, seeding = tr.seeding,
        units = tr.units, schedule = tr.schedule, days = tr.days[keep], mean = tr.mean[keep], sd = tr.sd[keep],
        sem = tr.sem[keep])
end

"Local dip at interior index i: log y_i minus the mean of its neighbours' logs (nearest remaining days)."
local_dip(tr, i) = log(max(tr.mean[i], 1e-3)) - (log(max(tr.mean[i - 1], 1e-3)) + log(max(tr.mean[i + 1], 1e-3))) / 2

const PLATE_DIP = -0.25     # plate-wide: median dip across the plate's series
const POINT_DEV = log(2.0)  # single series: interior point > 2-fold off its neighbours
const END_DEV = log(2.5)    # last point > 2.5-fold off the extrapolation of the previous two

"""
Plate-wide artefact days: on plates with ≥ 3 series, repeatedly remove the day whose median local
dip is most negative while it is below PLATE_DIP (recomputed after each removal, so two-day dips
such as LR run 1 days 8-9 are found).
"""
function qc_plate_days(trajs)
    for plate in unique(t.plate for t in trajs)
        idx = findall(t -> t.plate == plate, trajs)
        length(idx) >= 3 || continue
        while true
            dips = Dict{Float64, Vector{Float64}}()
            for k in idx, i in 2:length(trajs[k].days) - 1
                push!(get!(dips, trajs[k].days[i], Float64[]), local_dip(trajs[k], i))
            end
            isempty(dips) && break
            worst = argmin(d -> median(dips[d]), collect(keys(dips)))
            median(dips[worst]) < PLATE_DIP && length(dips[worst]) >= 3 || break
            for k in idx
                keep = trajs[k].days .!= worst
                all(keep) || (trajs[k] = keep_points(trajs[k], keep, "plate-wide artefact day ($(plate))"))
            end
        end
    end
    return trajs
end

_extrap(ly, days, i, j, k) = ly[j] + (ly[j] - ly[i]) * (days[k] - days[j]) / (days[j] - days[i])

"""
Last point > 2.5-fold off the extrapolation of the previous two, judged only when the point before
it is itself on trend (within 2-fold of its own extrapolation); otherwise that earlier point is the
anomaly and the interior rule handles it.
"""
function qc_last_point(tr::Trajectory)
    n = length(tr.days)
    n >= 5 || return tr
    ly = log.(max.(tr.mean, 1e-3)); d = tr.days
    prev_on_trend = abs(ly[n - 1] - _extrap(ly, d, n - 3, n - 2, n - 1)) < log(2.0)
    if prev_on_trend && abs(ly[n] - _extrap(ly, d, n - 2, n - 1, n)) > END_DEV
        keep = trues(n); keep[n] = false
        tr = keep_points(tr, keep, "last-point outlier (>2.5x off trend)")
    end
    return tr
end

"""
Single-series outliers, applied after the plate rule: last-point check, then isolated interior points
(> 2-fold above, or below, BOTH neighbours; at most two per series), then the last-point check again.
"""
function qc_series_points(tr::Trajectory)
    tr = qc_last_point(tr)
    for _ in 1:2
        length(tr.days) >= 4 || break
        ly = log.(max.(tr.mean, 1e-3))
        iso = [min(ly[i] - ly[i - 1], ly[i] - ly[i + 1]) > POINT_DEV || max(ly[i] - ly[i - 1], ly[i] - ly[i + 1]) < -POINT_DEV
               for i in 2:length(ly) - 1]
        any(iso) || break
        devs = [iso[i - 1] ? abs(local_dip(tr, i)) : 0.0 for i in 2:length(ly) - 1]
        j = argmax(devs)
        keep = trues(length(tr.days)); keep[j + 1] = false
        tr = keep_points(tr, keep, "isolated point (>2x vs both neighbours)")
    end
    return qc_last_point(tr)
end

# Day 8-9 rule (user, 2026-09-29: "only apply where drops exist"): in every series, day 8 and/or 9 is
# removed when it dips below the log-linear line joining the nearest remaining days before 8 and after 9.
const USER_DAYS = (8.0, 9.0)
const USER_DIP = -0.25

function qc_days_8_9(tr::Trajectory)
    before = findlast(<(8.0), tr.days); after = findfirst(>(9.0), tr.days)
    (before === nothing || after === nothing) && return tr
    l(i) = log(max(tr.mean[i], 1e-3))
    t0, t1 = tr.days[before], tr.days[after]
    keep = trues(length(tr.days))
    for (i, d) in enumerate(tr.days)
        d in USER_DAYS || continue
        line = l(before) + (l(after) - l(before)) * (d - t0) / (t1 - t0)
        l(i) - line < USER_DIP && (keep[i] = false)
    end
    return all(keep) ? tr : keep_points(tr, keep, "day 8-9 drop (below the day-7 to day-10 line)")
end

"Load every 20k and 30k trajectory, apply QC in order (plate days → single points → day 8-9 drops), check the hold-out."
function load_all()
    empty!(EXCLUDED_POINTS); empty!(EXCLUDED_WELLS)
    trajs = vcat(load_raw(), load_lr())
    trajs = qc_plate_days(trajs)
    trajs = [qc_series_points(t) for t in trajs]
    trajs = [qc_days_8_9(t) for t in trajs]
    for tr in trajs
        if tr.split != :train
            @assert tr.role == :validation "validation trajectory $(tr.id) leaked into a training role"
        end
        tr.split == :train && @assert tr.source != :ce1 "$(tr.id) must be held out"
        tr.split == :train && tr.role == :competition && @assert !treated(tr) "$(tr.id): treated data in competition role"
        @assert tr.seeding in (20, 30)
        @assert all(>=(1.0), tr.days)
    end
    @assert length(unique(t.id for t in trajs)) == length(trajs) "duplicate trajectory ids"
    return trajs
end

training(trajs; role = nothing) = filter(t -> t.split == :train && (role === nothing || t.role == role), trajs)
validation(trajs) = filter(t -> t.split != :train, trajs)

"Well-level co-culture pairs (S, R) keyed by group."
function coculture_pairs(trajs)
    pairs = Dict{String, NamedTuple}()
    for g in unique(t.group for t in trajs if t.culture == :co)
        members = filter(t -> t.group == g, trajs)
        s = findfirst(t -> t.line == :S, members)
        r = findfirst(t -> t.line == :R, members)
        (s === nothing || r === nothing) && continue
        pairs[g] = (S = members[s], R = members[r])
    end
    return pairs
end

function long_table(trajs)
    rows = DataFrame()
    for tr in trajs, (i, d) in enumerate(tr.days)
        push!(rows, (; traj_id = tr.id, source = String(tr.source), plate = tr.plate, split = String(tr.split),
            role = String(tr.role), line = String(tr.line), resource = String(tr.resource),
            culture = String(tr.culture), group = tr.group, s_frac = tr.s_frac, seeding_k = tr.seeding,
            units = String(tr.units), schedule = describe(tr.schedule), dose = dose_at(tr.schedule, d),
            day = d, mean = tr.mean[i], sd = tr.sd[i], sem = tr.sem[i]))
    end
    return rows
end

"Filtered drug exposure v(t) with dv/dt = gamma (dose(t) - v), v(0) = 0 (exact, piecewise)."
function exposure(s::Schedule, gamma, t)
    v, tprev = zero(gamma * t), 0.0
    for (i, st) in enumerate(s.starts)
        st >= t && break
        tend = i < length(s.starts) ? min(s.starts[i + 1], t) : t
        v = s.doses[i] + (v - s.doses[i]) * exp(-gamma * (tend - max(st, tprev)))
        tprev = tend
    end
    return v
end

"Cumulative dose (uM * day) up to t."
cumdose(s::Schedule, t) = sum(d * max(0.0, min(t, i < length(s.starts) ? s.starts[i + 1] : Inf) - st)
                              for (i, (st, d)) in enumerate(zip(s.starts, s.doses)))
since_switch(s::Schedule, t) = t - s.starts[max(searchsortedlast(s.starts, t), 1)]

"""
Global area -> count factor for intermittent data: ratio of intermittent untreated-first
(off-on) day-7 value to untreated monoculture day-7 value, per line, averaged.
Only used for the next-day feature scale; the ODE models profile a per-trajectory offset instead.
"""
function area_per_count(trajs)
    ratios = Float64[]
    for line in (:S, :R)
        a = filter(t -> t.source == :intermittent && t.line == line && t.seeding == 30 && t.schedule.doses[1] == 0.0, trajs)
        c = filter(t -> t.source == :mono && t.line == line && t.seeding == 30 && !treated(t), trajs)
        (length(a) == 1 && length(c) == 1) || continue
        i, j = findfirst(==(7.0), a[1].days), findfirst(==(7.0), c[1].days)
        (i === nothing || j === nothing) && continue
        push!(ratios, a[1].mean[i] / c[1].mean[j])
    end
    return isempty(ratios) ? 1.0 : mean(ratios)
end

"Observed values on the count scale."
count_scale(tr::Trajectory, apc) = tr.units == :area ? tr.mean ./ apc : tr.mean
