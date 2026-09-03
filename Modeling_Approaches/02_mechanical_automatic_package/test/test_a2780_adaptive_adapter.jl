using .MechanicalAutomaticModeling.A2780AdaptiveAdapter
using CSV
using DataFrames
using GrowthParameterEstimation

@testset "A2780 fitted-artifact adapter" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    config = load_a2780_adaptive_config(package_root)
    stage3_status = CSV.read(joinpath(package_root, "outputs", "csv", "coculture_untreated", "coculture_untreated_pooling_status.csv"), DataFrame)
    stage4_status = CSV.read(joinpath(package_root, "outputs", "csv", "coculture_treated", "linked_treatment_status.csv"), DataFrame)
    @test config.winners.stage3.model == String(first(stage3_status.winning_model))
    @test config.winners.stage4.model == String(first(stage4_status.winning_model))
    @test config.winners.stage4.eligible
    @test config.observables == ["A2780Naive", "Total A2780cis", "Total population"]
    @test "cis tolerant-like" in config.latent_states_hidden
    @test all(isfinite, values(config.winners.stage4.parameters))
    @test any(candidate -> candidate.stage == 1 && candidate.cell_line == "A2780Naive" && candidate.eligible, config.advanced_models)
    @test any(candidate -> candidate.stage == 1 && candidate.cell_line == "A2780cis" && candidate.eligible, config.advanced_models)
    @test any(candidate -> candidate.stage == 2 && candidate.eligible, config.advanced_models)
    @test any(candidate -> candidate.stage == 3 && candidate.eligible, config.advanced_models)
    @test any(candidate -> candidate.stage == 4 && candidate.eligible, config.advanced_models)
    @test length(config.experimental_library.original) == 6
    @test length(config.experimental_library.fitted) == 6
    @test haskey(config.experimental_library.original, "coculture_untreated__20k__50-50__0")
    @test haskey(config.experimental_library.original, "coculture_untreated__30k__50-50__0")
    untreated_20k = config.experimental_library.fitted["coculture_untreated__20k__50-50__0"]
    untreated_30k = config.experimental_library.fitted["coculture_untreated__30k__50-50__0"]
    @test untreated_20k.naive[1].t == 0
    @test untreated_20k.naive[1].y + untreated_20k.cis[1].y == 67
    @test untreated_30k.naive[1].y + untreated_30k.cis[1].y == 100

    zero_protocol = TreatmentProtocol([TreatmentWindow(0, 14, 0, 2)]; horizon = 14)
    zero_scenario = build_a2780_scenario(config, zero_protocol; initial_state = [50.0, 25.0, 25.0, 0.0, 0.0])
    @test zero_scenario.model isa GrowthParameterEstimation.ModelSpec
    zero_result = simulate_protocol(zero_scenario)
    @test all(state -> all(>=(0), state), zero_result.states)
    @test all(==(0), zero_result.effective_exposure)
    @test length(zero_result.observables[1]) == 3
    @test zero_result.observables[1] == [50.0, 50.0, 100.0]
    zero_summary = summarize_outcomes(zero_result, (_, observable) -> observable[3] > 1e9)
    @test zero_summary.median_final_total ≈ zero_result.observables[end][3]

    overlay = CSV.read(joinpath(package_root, "outputs", "csv", "coculture_untreated", "figures", "coculture_untreated_joint_overlays.csv"), DataFrame)
    reference = filter(row ->
        row.model == "lv_asymmetric_competition_death" &&
        row.pooling_mode == "partial_5pct" &&
        row.density == "20k" &&
        row.mix == "50-50",
        overlay,
    )
    exported_naive = sort(filter(row -> row.component == "sensitive", reference), :time)
    exported_cis = sort(filter(row -> row.component == "resistant", reference), :time)
    reproduction = simulate_protocol(build_a2780_scenario(config, zero_protocol;
        density = "20k",
        initial_state = [33.5, 16.75, 16.75, 0.0, 0.0],
    ))
    for row in eachrow(exported_naive)
        index = findfirst(time -> isapprox(time, row.time; atol = 1e-8), reproduction.times)
        @test isapprox(reproduction.observables[index][1], row.predicted; rtol = 0.01)
    end
    for row in eachrow(exported_cis)
        index = findfirst(time -> isapprox(time, row.time; atol = 1e-8), reproduction.times)
        @test isapprox(reproduction.observables[index][2], row.predicted; rtol = 0.01)
    end

    linked_overlay = CSV.read(joinpath(package_root, "outputs", "csv", "coculture_treated", "figures", "linked_treatment_combined_overlays.csv"), DataFrame)
    linked_reference = filter(row ->
        row.model == config.winners.stage4.model &&
        row.context == "coculture" &&
        row.density == "20k" &&
        row.mix == "50-50",
        linked_overlay,
    )
    fitted_fraction = build_a2780_scenario(config, zero_protocol).parameters.initial_tolerant_fraction
    fitted_initial_state = [33.5, 33.5 * (1 - fitted_fraction), 33.5 * fitted_fraction, 0.0, 0.0, 0.0]
    fitted_protocol = TreatmentProtocol([TreatmentWindow(0, 14, 1.0, 2)]; horizon = 14)
    fitted_treated = simulate_protocol(build_a2780_scenario(config, fitted_protocol;
        density = "20k",
        initial_state = fitted_initial_state,
    ))
    for component in ("sensitive", "resistant")
        observable_index = component == "sensitive" ? 1 : 2
        for row in eachrow(sort(filter(row -> row.component == component, linked_reference), :time))
            index = findfirst(time -> isapprox(time, row.time; atol = 1e-8), fitted_treated.times)
            @test isapprox(fitted_treated.observables[index][observable_index], row.predicted; rtol = 0.01)
        end
    end

    treated_protocol = TreatmentProtocol([TreatmentWindow(0, 4, 1.0, 2)]; horizon = 8)
    treated = simulate_protocol(build_a2780_scenario(config, treated_protocol))
    @test maximum(treated.effective_exposure) > 0
    @test treated.effective_exposure[end] < maximum(treated.effective_exposure)
end
