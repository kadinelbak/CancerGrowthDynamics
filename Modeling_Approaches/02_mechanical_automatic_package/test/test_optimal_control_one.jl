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
end
