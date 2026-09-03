using DataFrames
using CSV

@testset "Optimal Control One" begin
    OC = MechanicalAutomaticModeling.OptimalControlOne
    @test OC.A2780_DOSE_MAP["IC25"] == 1.47
    @test OC.A2780_DOSE_MAP["IC50"] == 1.0
    @test OC.A2780_DOSE_MAP["IC75"] == 0.67
    p = OC.PaperParameters(0.4, 0.7, 0.2, 0.15, 0.1)
    times = collect(0.0:0.5:4.0)
    states = OC.simulate_two_population(p, times)
    @test size(states) == (length(times), 2)
    @test states[1, :] == [1.0, 0.0]
    @test all(states .>= 0)
    @test states[end, 2] > 0

    three = OC.simulate_three_population([0.4, 0.7, 0.2, 0.3, 0.15, 0.1], times)
    @test size(three) == (length(times), 3)
    @test all(three .>= 0)
    @test three[end, 2] > 0
    @test three[end, 3] > 0

    fitted = Dict(0.67 => p, 1.47 => OC.PaperParameters(0.6, 0.9, 0.1, 0.12, 0.08))
    middle = OC.interpolate_parameters(1.0, fitted)
    @test p.rS < middle.rS < 0.6
    @test p.alpha > middle.alpha > 0.1

    doses = [0.67, 1.0, 1.47]
    for name in OC.CANDIDATE_MODELS
        definition = OC._candidate_definition(name, length(doses))
        candidate = [(lo + hi) / 2 for (lo, hi) in definition.bounds]
        prediction = OC._simulate_candidate(name, candidate, doses, 1.0, times)
        @test length(candidate) == length(definition.names)
        @test length(prediction.total) == length(times)
        @test all(isfinite, prediction.total)
        @test all(prediction.total .>= 0)
    end

    selected_definition = OC._candidate_definition("shared-growth delayed logistic", 3)
    selected_fit = (name = "shared-growth delayed logistic", params = [0.5, 8.0, 0.1, 0.3, 0.8], doses = doses)
    control = OC._simulate_selected_control(selected_fit, fill(1.0, 4); horizon = 4.0)
    @test all(isfinite, control.total)
    @test all(control.total .>= 0)
    @test control.latent == zeros(length(control.time))

    ranking = DataFrame(
        model = ["simple", "complex"], pooling_mode = ["shared", "shared"],
        n_parameters = [2, 4], bic = [10.0, 12.0], boundary_issue = [false, true],
        eligible_for_inheritance = [true, true],
    )
    status = DataFrame(winning_model = ["simple"], winning_pooling_mode = ["shared"])
    candidates = OC._top_stage_rows(ranking, status)
    @test candidates.candidate == ["M1", "M2"]
    @test candidates.selected == [true, false]

    overlay = DataFrame(
        model = repeat(["simple", "complex"], inner = 3), pooling_mode = fill("shared", 6),
        observed = repeat([1.0, 2.0, 3.0], 2), predicted = [1.0, 2.0, 3.0, 1.0, 1.0, 1.0],
        density = fill("30k", 6),
    )
    audited = OC._attach_candidate_errors(candidates, overlay)
    @test audited.median_nRMSE[1] == 0.0
    @test audited.median_nRMSE[2] > audited.median_nRMSE[1]

    report = read(normpath(joinpath(@__DIR__, "..", "outputs", "reports", "optimal_control_one.html")), String)
    @test contains(report, "Conditional Model Tournament")
    @test findfirst("Conditional Model Tournament", report) < findfirst("Stage 1: Untreated Monoculture", report)
    @test all(contains(report, "Stage $stage:") for stage in 1:4)
    @test contains(report, "The optimizer was run despite incomplete validation")
    @test contains(report, "Why The Parameters Are Not Yet Reliable")
    @test contains(report, "Five rows pairing each 24-day cisplatin schedule")
    @test contains(report, "Days 14-24 are forward model projections")
    @test contains(report, "Stage 4 Candidate Endpoint Audit")
    @test !contains(report, "endpoint-only prediction")

    control = CSV.read(normpath(joinpath(@__DIR__, "..", "outputs", "csv", "optimal_control_one", "control_metrics.csv")), DataFrame)
    @test Set(control.schedule) == Set(["Optimized", "Constant 1 uM", "Front-loaded", "Back-loaded", "Pulsed"])
    @test maximum(control.dose_AUC) - minimum(control.dose_AUC) < 1e-8
    @test isapprox(only(filter(:schedule => ==("Constant 1 uM"), control).dose_AUC), 24.0; atol = 1e-8)
    @test only(filter(:schedule => ==("Optimized"), control).final_total) <= only(filter(:schedule => ==("Constant 1 uM"), control).final_total)

    trajectories = CSV.read(normpath(joinpath(@__DIR__, "..", "outputs", "csv", "optimal_control_one", "control_trajectories.csv")), DataFrame)
    @test isapprox(maximum(trajectories.time_day), 24.0; atol = 1e-8)
end
