# Stage 0: build, quality-control, audit and plot every 20k and 30k trajectory.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots
const N = NNTrack

trajs = N.load_all()
out = N.outdir("data")
tbl = N.long_table(trajs)
CSV.write(joinpath(out, "harmonised_long.csv"), tbl)
CSV.write(joinpath(out, "excluded_points.csv"), N.EXCLUDED_POINTS)
CSV.write(joinpath(out, "excluded_well_days.csv"), N.EXCLUDED_WELLS)

summary = combine(groupby(tbl, [:split, :role, :source, :plate, :line, :resource, :culture]),
    :traj_id => (x -> length(unique(x))) => :trajectories, nrow => :points)
sort!(summary, [:split, :role, :source, :plate, :line])
CSV.write(joinpath(out, "summary.csv"), summary)
println(summary)
rules = sort(combine(groupby(N.EXCLUDED_POINTS, :rule), nrow => :points, :traj_id => (x -> length(unique(x))) => :trajectories), :points; rev = true)
println(rules)

# ---- checks ----
@assert count(t -> t.source == :mono, trajs) == 16                  # 2 densities × 2 lines × (0, 0.62, 1, 1.47 uM)
@assert count(t -> t.source == :intermittent, trajs) == 8           # 2 densities × 2 lines × 2 schedules
@assert count(t -> t.source == :ce0, trajs) == 12 && count(t -> t.source == :ce1, trajs) == 12
@assert count(t -> t.source == :lr, trajs) == 22                    # run 1: 2 densities × (3 ratios × 2 lines + 2 mono); run 2: 30k, 3 × 2
@assert all(t.seeding in (20, 30) for t in trajs)

# ---- inventory of every data file considered (audit 2026-09-28) ----
inventory = DataFrame(
    data = ["Untreated monoculture 20k + 30k (Processed_Datasets/Untreated MonoCulture/{20k,30k})",
            "Treated monoculture 20k + 30k, IC25/IC50/IC75 folders = 1.47 / 1.0 / 0.62 uM",
            "Intermittent 20k + 30k (on→off, off→on at 1 uM)",
            "Untreated co-culture, raw flat files (legacy), plate columns 1-3 = 20k, 4-6 = 30k",
            "Treated co-culture (1 uM), raw flat files (legacy), plate columns 1-3 = 20k, 4-6 = 30k",
            "LR workbook: run 1 (20k + 30k), run 2 (30k)",
            "New Datasets A2780*_30k.csv (8 files)",
            "New Datasets ce0_*/ce1_* workbooks",
            "Processed_Datasets co-culture 20k/30k Averages",
            "TykNu (all densities); 40k",
            "New Datasets/Validation (LR + cisplatin)"],
    status = ["USED (rebuilt from per-image data)", "USED (rebuilt from per-image data)", "USED (rebuilt from per-image data)",
              "USED (rebuilt by well column)", "USED as held-out validation (rebuilt by well column)", "USED (no raw images exist)",
              "not loaded: identical to the raw monoculture day means (max rel. diff < 1e-6); day 0 = 100 is a normalisation",
              "not loaded: identical to the column-4-6 rebuild, except ce0 75S25R resistant (~20% lower, different threshold) and one day-6 value",
              "not loaded: split by imaging tile instead of well column, so each 'density' mixes 20k and 30k wells; the same raw files are re-split correctly by column above",
              "excluded (out of scope)",
              "excluded (user decision)"],
    note = ["A2780cis runs to day 18 (days 16, 18 are extra vs New Datasets); wells: row A = 20k, B = 30k", "3 wells × 6 tiles per day; rows A/C = 20k, B/D = 30k", "rows A/C = 20k, B/D = 30k", "3 wells per ratio and density", "3 wells per ratio and density",
            "day 0 is the seeding value, dropped", "", "", "", "", ""])
CSV.write(joinpath(out, "inventory.csv"), inventory)

# ---- plots: one panel per trajectory, same cell-count axis, shared key; excluded points by rule ----
gr()
rule_style = Dict("plate" => (:red, "plate-wide artefact day"), "user" => (:orange, "day 8-9 drop"),
                  "isolated" => (:purple, "isolated point"), "last" => (:magenta, "last-point outlier"))
rulekey(r) = startswith(r, "plate") ? "plate" : startswith(r, "day 8-9") ? "user" : startswith(r, "isolated") ? "isolated" : "last"
function data_grid(sel, fname; ncol = 4)
    sub = filter(sel, trajs)
    ex = N.EXCLUDED_POINTS[in.(N.EXCLUDED_POINTS.traj_id, Ref(Set(t.id for t in sub))), :]
    allv = vcat(vcat((t.mean for t in sub)...), ex.mean)
    yl = (10.0^floor(log10(max(minimum(allv), 1.0))), 10.0^ceil(log10(maximum(allv))))
    ps = []
    for t in sub
        p = plot(; title = "$(t.id)\n$(N.describe(t.schedule))", titlefontsize = 8, legend = false, yscale = :log10, ylims = yl,
                 xlims = (0, 19), xlabel = "day", ylabel = "cells", guidefontsize = 7, tickfontsize = 6)
        length(t.schedule.doses) > 1 && vline!(p, [7.0]; color = :gray, ls = :dash)
        plot!(p, t.days, t.mean; yerror = nothing, color = t.line == :S ? :crimson : :royalblue, lw = 1.5, marker = :circle, ms = 3)
        e = ex[ex.traj_id .== t.id, :]
        for r in eachrow(e)
            c, _ = rule_style[rulekey(r.rule)]
            scatter!(p, [r.day], [max(r.mean, 1.0)]; marker = :circle, ms = 5, mc = :white, msc = c, msw = 2)
        end
        push!(ps, p)
    end
    key = plot(; framestyle = :none, legend = :top, legend_column = 3, legendfontsize = 11, foreground_color_legend = nothing, xlims = (0, 1), ylims = (0, 1))
    plot!(key, [-1, -1], [-1, -1]; color = :crimson, lw = 4, label = "A2780 (S)")
    plot!(key, [-1, -1], [-1, -1]; color = :royalblue, lw = 4, label = "A2780cis (R)")
    for k in ("plate", "user", "isolated", "last")
        c, lab = rule_style[k]
        scatter!(key, [-1], [-1]; mc = :white, msc = c, msw = 2, ms = 7, label = "excluded: $lab")
    end
    nrow_ = cld(length(ps), ncol)
    savefig(plot(key, plot(ps...; layout = (nrow_, ncol)); layout = @layout([a{0.08h}; b]), size = (1500, 150 + 300 * nrow_),
                 left_margin = 5Plots.mm, bottom_margin = 3Plots.mm), joinpath(out, fname))
    return fname
end
figs = [("Monoculture (untreated and constant dose)", data_grid(t -> t.source == :mono, "traj_mono.png")),
        ("Intermittent dosing (dashed line = switch at day 7)", data_grid(t -> t.source == :intermittent, "traj_intermittent.png")),
        ("Untreated co-culture (ce0)", data_grid(t -> t.source == :ce0, "traj_ce0.png"; ncol = 6)),
        ("Treated co-culture (ce1, held out)", data_grid(t -> t.source == :ce1, "traj_ce1.png"; ncol = 6)),
        ("Low resource runs 1 and 2", data_grid(t -> t.source == :lr, "traj_lr.png"))]
# density comparison: 20k vs 30k for the same condition
dens = combine(groupby(filter(r -> r.day in (1.0, 7.0, 14.0), tbl), [:source, :line, :schedule, :s_frac, :resource, :day]),
    [:seeding_k, :mean] => ((k, m) -> (any(k .== 20) && any(k .== 30)) ? m[findfirst(==(30), k)] / m[findfirst(==(20), k)] : NaN) => :ratio_30k_over_20k)
dens = dropmissing(filter(r -> isfinite(r.ratio_30k_over_20k), dens))
CSV.write(joinpath(out, "density_ratio.csv"), dens)

N.write_report(joinpath(out, "report.html"), "Stage 0 — data audit and harmonised data",
    [("Inventory: every data source in the repo", N.para("All 20k and 30k A2780 / A2780cis data in the repo is used. Where raw per-image files exist, day means are rebuilt from per-well values so bad wells can be caught. Pre-averaged files that duplicate the raw data are not loaded twice.") * N.html_table(inventory)),
     ("Quality control, in order", N.para("""
        <b>1. Well level (raw sources):</b> a well-day is dropped if that well jumps > 2-fold against its own previous/next day while the other wells move < 1.5-fold (a single bad image set). Wells that are consistently high or low (seeding differences) are kept. $(nrow(N.EXCLUDED_WELLS)) well-days removed; see <code>excluded_well_days.csv</code>.<br>
        <b>2. Plate-wide artefact days:</b> on each plate/experiment with ≥ 3 series, a day is removed for every series when the median dip against neighbouring days is below $(N.PLATE_DIP) (log), re-checked after each removal so two-day dips are found.<br>
        <b>3. Single points:</b> a last point > 2.5-fold off trend (only if the point before is on trend), and interior points > 2-fold above or below both neighbours (at most two per series).<br>
        <b>4. Day 8-9 drops:</b> in every series, day 8 and/or 9 is removed only when it lies more than $(-N.USER_DIP) (log) below the line joining the nearest remaining days before day 8 and after day 9.""") *
        N.html_table(rules) * N.html_table(N.EXCLUDED_WELLS)),
     ("Splits", N.para("ce1 treated co-culture is held out; the loader asserts it never appears in training. Untreated co-culture is used only for competition coefficients.") * N.html_table(summary)),
     ("20k vs 30k", N.para("Ratio of 30k to 20k day means for the same condition at days 1, 7 and 14. Seeding ratio is 1.5; a ratio falling towards 1 by day 14 means both densities approach the same carrying capacity (the model's initial condition carries density; K and rates are shared unless the fits say otherwise).") * N.html_table(dens)),
     ("Trajectories (open circles = excluded points, coloured by rule)", join(N.figure(f, cap) for (cap, f) in figs)),
     ("Assumptions", N.para("Intermittent: <code>T</code> = 1 uM days 1–7 then off; <code>UT</code> = off then 1 uM days 8–14 (Understanding.md, Well Layouts). Treated monoculture folders IC25 / IC50 / IC75 = 1.47 / 1.0 / 0.62 uM (identical to the New Datasets files named with those doses). Co-culture ratios are sensitive:resistant at seeding; column 4 = 75/25, 5 = 50/50, 6 = 25/75."))])
println("wrote ", joinpath(out, "report.html"))
