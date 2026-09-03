using .MechanicalAutomaticModeling.AdaptiveSimulationEngine

@testset "AdaptiveSimulationEngine" begin
    holiday = TreatmentWindow(2, 4, 0, 2; label = "holiday")
    first = TreatmentWindow(0, 2, 1, 1; label = "first")
    second = TreatmentWindow(2, 4, 0.5, 1; label = "second")
    @test TreatmentProtocol([second, first, holiday]; horizon = 5).windows[1] == first
    @test_throws ArgumentError TreatmentProtocol([
        TreatmentWindow(0, 3, 1, 1),
        TreatmentWindow(2, 4, 1, 1),
    ])
    @test_throws ArgumentError TreatmentProtocol([
        TreatmentWindow(0, 4, 1, 1),
        TreatmentWindow(1, 2, 0, 1; label = "holiday"),
        TreatmentWindow(3, 5, 1, 1),
    ])

    callback = (_, state, parameters, t0, t1, dose) -> (
        state = max.(0.0, state .+ (t1 - t0) .* (parameters.growth .* state .- dose .* parameters.kill .* state)),
        effective_exposure = dose,
    )
    protocol = TreatmentProtocol([first, second, holiday]; monitoring_interval = 1, horizon = 5)
    scenario = SimulationScenario(:linear, (growth = [0.1, 0.05], kill = [0.2, 0.1]), [100.0, 50.0], protocol;
        observable = identity,
        segment_simulator = callback,
        save_interval = 0.25,
    )
    result = simulate_protocol(scenario)
    @test result.times[1] == 0
    @test result.times[end] == 5
    @test all(state -> all(>=(0), state), result.states)
    @test any(==(1.0), result.commanded_dose)
    @test any(==(0.0), result.commanded_dose)

    decision_times = Float64[]
    policy = context -> begin
        push!(decision_times, context.time)
        context.time < 2 ? 1.0 : 0.0
    end
    adaptive = TreatmentProtocol(TreatmentWindow[];
        monitoring_interval = 1,
        decision_policy = policy,
        allowable_doses = [0.0, 1.0],
        horizon = 3,
    )
    adaptive_result = simulate_protocol(SimulationScenario(:linear, scenario.parameters, [100.0, 50.0], adaptive;
        observable = identity,
        segment_simulator = callback,
        save_interval = 0.2,
    ))
    @test decision_times == [0.0, 1.0, 2.0]
    @test adaptive_result.commanded_dose[end] == 0

    ensemble1 = simulate_ensemble(scenario, [scenario.parameters, scenario.parameters])
    ensemble2 = simulate_ensemble(scenario, [scenario.parameters, scenario.parameters])
    @test ensemble1[2].states == ensemble2[2].states
    summary = summarize_outcomes(ensemble1, (_, observable) -> sum(observable) > 10_000)
    @test summary.probability_of_control == 1
    @test summary.median_final_total > 0
    @test summary.members[1].cumulative_dose == 3.0
    @test summary.members[1].cycle_count == 1

    reference = [(time = t, state = [100.0 + 10t, 50.0]) for t in 0:5]
    lookup = equivalent_untreated_time([131.0, 50.0], reference)
    @test lookup.time == 3
    @test equivalent_untreated_time([1.0], NamedTuple[]) === nothing
end
