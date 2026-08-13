using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using Test

root = @__DIR__
csv_root = joinpath(root, "outputs", "csv")
image_root = joinpath(root, "outputs", "images")
report_path = joinpath(root, "outputs", "reports", "a2780_staged_model_comparison.html")

finite_metric(values) = all(value -> isfinite(Float64(value)) && abs(Float64(value)) < 9.99e11, values)

@testset "A2780 timing and linked-treatment validation" begin
    timing_path = joinpath(csv_root, "monoculture_treated", "monoculture_treated_timing_hypothesis_ranking.csv")
    timing = sort!(CSV.read(timing_path, DataFrame), :bic)
    expected_timing = Set([
        "independent_onset_gradual",
        "shared_onset_gradual",
        "partial_onset_0_5day",
        "cis_gradual_only",
        "cis_onset_only",
    ])
    @test nrow(timing) == 5
    @test Set(String.(timing.model)) == expected_timing
    @test finite_metric(timing.bic)
    @test finite_metric(timing.scaled_ssr)
    @test String(first(timing.model)) == "cis_onset_only"
    @test all(String.(timing.cis_activation_mode[String.(timing.model) .== "cis_onset_only"]) .== "step")

    timing_starts = CSV.read(joinpath(csv_root, "monoculture_treated", "monoculture_treated_timing_hypothesis_multistart.csv"), DataFrame)
    @test nrow(timing_starts) == 15
    @test all(String.(timing_starts.status) .== "completed")
    @test finite_metric(timing_starts.bic)

    timing_overlay = CSV.read(joinpath(csv_root, "monoculture_treated", "figures", "monoculture_treated_timing_hypothesis_overlays.csv"), DataFrame)
    timing_panels = unique(select(timing_overlay, [:timing_hypothesis, :cell_line, :density, :dose]))
    @test nrow(timing_panels) == 60
    winner_timing_overlay = timing_overlay[String.(timing_overlay.timing_hypothesis) .== "cis_onset_only", :]
    @test nrow(unique(select(winner_timing_overlay, [:cell_line, :density, :dose]))) == 12
    for row in eachrow(winner_timing_overlay[iszero.(Float64.(winner_timing_overlay.time)), :])
        @test isapprox(Float64(row.observed), String(row.density) == "20k" ? 67.0 : 100.0; atol = 1e-8)
        @test isapprox(Float64(row.predicted), Float64(row.observed); atol = 1e-8)
    end

    linked = sort!(CSV.read(joinpath(csv_root, "coculture_treated", "linked_treatment_model_ranking.csv"), DataFrame), :bic)
    @test nrow(linked) == 9
    @test finite_metric(linked.bic)
    @test finite_metric(linked.scaled_ssr)
    @test all(String.(linked.inherited_timing_hypothesis) .== "cis_onset_only")
    @test String(first(linked.model)) == "load_plus_tolerant_growth_context"
    @test Bool(first(linked.eligible_for_inheritance))

    status = first(CSV.read(joinpath(csv_root, "coculture_treated", "linked_treatment_status.csv"), DataFrame))
    @test String(status.winning_model) == "load_plus_tolerant_growth_context"
    @test Bool(status.inheritance_allowed)
    @test !Bool(status.inadequate_inheritance)
    @test Float64(status.fully_free_bic_improvement) < 10.0

    audit = CSV.read(joinpath(csv_root, "coculture_treated", "linked_treatment_shared_parameter_audit.csv"), DataFrame)
    inherited = audit[isfinite.(Float64.(audit.sequential_stage2_seed)), :]
    @test maximum(abs.(Float64.(inherited.deviation_percent))) <= 5.000001
    growth_row = only(eachrow(audit[String.(audit.parameter) .== "log_cis_tolerant_growth_context", :]))
    @test 1.0 < exp(Float64(growth_row.estimate)) < 2.0

    linked_overlay = CSV.read(joinpath(csv_root, "coculture_treated", "figures", "linked_treatment_combined_overlays.csv"), DataFrame)
    winner_overlay = linked_overlay[String.(linked_overlay.model) .== String(status.winning_model), :]
    mono = winner_overlay[String.(winner_overlay.context) .== "monoculture", :]
    co = winner_overlay[String.(winner_overlay.context) .== "coculture", :]
    @test nrow(unique(select(mono, [:cell_line, :density, :dose]))) == 12
    @test nrow(unique(select(co, [:density, :mix, :component]))) == 12
    for group in groupby(co[Bool.(co.fixed_day0_anchor), :], [:density, :mix, :model])
        expected = String(first(group.density)) == "20k" ? 67.0 : 100.0
        @test isapprox(sum(Float64.(group.observed)), expected; atol = 1e-8)
        @test isapprox(sum(Float64.(group.predicted)), expected; atol = 1e-8)
    end

    overview = CSV.read(joinpath(root, "outputs", "reports", "a2780_stage_overview.csv"), DataFrame)
    @test String.(overview.condition) == [
        "monoculture_untreated",
        "monoculture_treated",
        "coculture_untreated",
        "coculture_treated",
    ]
    @test all(String.(overview.status) .== "completed")

    source = read(joinpath(root, "src", "timing_hypotheses.jl"), String) *
        read(joinpath(root, "src", "linked_treatment_joint.jl"), String)
    @test occursin("GrowthParameterEstimation.run_joint_multistart", source)
    @test occursin("GrowthParameterEstimation.run_joint_fit", source)
    @test occursin("profile_joint_fit_bounds_two_sided", source)

    for path in [
        joinpath(image_root, "monoculture_treated", "figures", "monoculture_treated_timing_hypothesis_grid.png"),
        joinpath(image_root, "coculture_treated", "figures", "linked_treatment_monoculture_grid.png"),
        joinpath(image_root, "coculture_treated", "figures", "linked_treatment_coculture_grid.png"),
    ]
        @test isfile(path)
        @test filesize(path) > 50_000
    end

    html = read(report_path, String)
    @test occursin("Onset and gradual-effect combinations", html)
    @test occursin("load_plus_tolerant_growth_context", html)
    @test occursin("tex-chtml.js", html)
    @test count("notation key", html) == 4
    @test occursin(raw"\frac{dX}{dt}", html)
    @test occursin(raw"\frac{dP}{dt}", html)
    @test occursin(raw"\frac{dN}{dt}", html)
    @test !occursin(raw"\frac{dL}{dt}", html)
    @test occursin("Different letters mean distinct biological states", html)
    @test occursin(raw"f_{T0}^{\mathrm{co}}", html)
    @test occursin(raw"\rho_T", html)
    @test occursin(raw"\mathrm{BIC}=n\ln", html)
    @test !occursin("<code>dN/dt", html)
    @test occursin("<table", html)
    @test count("<img", html) >= 5
end
