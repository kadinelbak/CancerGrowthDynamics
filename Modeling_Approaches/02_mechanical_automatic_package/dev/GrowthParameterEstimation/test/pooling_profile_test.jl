@testset "Pooling summaries and boundary profiles" begin
    pair = symmetric_relative_pair(10.0, log(1.05))
    @test isapprox(pair.low, 10 / 1.05)
    @test isapprox(pair.high, 10 * 1.05)

    pooling = DataFrame(
        cell_line = ["A", "A", "A", "B", "B", "B"],
        model = ["m", "m", "m", "m", "m", "m"],
        pooling_mode = ["shared", "partial_5pct", "independent_diagnostic", "shared", "partial_5pct", "independent_diagnostic"],
        bic = [100.0, 101.0, 95.0, 120.0, 118.0, 105.0],
        eligible_for_inheritance = [true, true, false, true, true, false],
    )
    summary = summarize_pooling_bic(pooling; top_n = 2)
    @test nrow(summary.ranking) == 4
    @test !summary.status[summary.status.cell_line .== "A", :inadequate_pooling][1]
    @test summary.status[summary.status.cell_line .== "B", :inadequate_pooling][1]

    function pooled_logistic!(du, u, p, t)
        r, K = p
        for i in eachindex(u)
            du[i] = r * u[i] * (1 - u[i] / K)
        end
    end
    times = collect(0.0:1.0:8.0)
    truth = [0.3, 25.0]
    initial = [1.0, 2.0]
    problem = ODEProblem(pooled_logistic!, initial, (first(times), last(times)), truth)
    solution = solve(problem, Tsit5(); saveat = times)
    datasets = [
        (x = times, y = [u[1] for u in solution.u], state_index = 1, residual_scale = 10.0),
        (x = times, y = [u[2] for u in solution.u], state_index = 2, residual_scale = 15.0),
    ]
    profiled = profile_joint_fit_bounds(
        pooled_logistic!, datasets, initial, [0.15, 20.0];
        bounds = [(0.01, 0.2), (5.0, 60.0)],
        parameter_names = [:r, :K],
        explicit_upper_profiles = Dict(:r => [0.3, 0.4]),
        optimizer = :nelder_mead,
        maxiters = 1200,
        show_stats = false,
    )
    @test isfinite(profiled.fit.bic)
    @test any(profiled.profile.accepted)
    @test profiled.bounds[1][2] > 0.2
    @test all(in.(
        profiled.identifiability.identifiability,
        Ref(["interior", "poorly_identified_at_bound"]),
    ))

    function lower_profile_model!(du, u, p, t)
        du[1] = p[1] * u[1]
    end
    lower_truth = [-0.5]
    lower_problem = ODEProblem(lower_profile_model!, [5.0], (0.0, 5.0), lower_truth)
    lower_solution = solve(lower_problem, Tsit5(); saveat = 0.0:1.0:5.0)
    lower_datasets = [(x = collect(0.0:1.0:5.0), y = [state[1] for state in lower_solution.u], state_index = 1, residual_scale = 5.0)]
    lower_profiled = profile_joint_fit_bounds_two_sided(
        lower_profile_model!, lower_datasets, [5.0], [-0.1];
        bounds = [(-0.2, 0.5)],
        parameter_names = [:beta],
        explicit_lower_profiles = Dict(:beta => [-0.6, -1.0]),
        physical_lower_limits = Dict(:beta => -1.0),
        optimizer = :nelder_mead,
        maxiters = 800,
        show_stats = false,
    )
    @test isfinite(lower_profiled.fit.bic)
    @test any((lower_profiled.profile.direction .== "lower") .& lower_profiled.profile.accepted)
    @test lower_profiled.bounds[1][1] < -0.2
end

@testset "Synthetic density pooling behavior" begin
    times = collect(0.0:1.0:6.0)
    base_u0 = [1.0, 1.5]

    function generate_pair(r, log_contrast; u0 = base_u0)
        function truth_model!(du, u, p, t)
            center, contrast = p
            du[1] = center * exp(-contrast) * u[1]
            du[2] = center * exp(contrast) * u[2]
        end
        problem = ODEProblem(truth_model!, u0, (first(times), last(times)), [r, log_contrast])
        solution = solve(problem, Tsit5(); saveat = times, reltol = 1e-10, abstol = 1e-10)
        return [[state[i] for state in solution.u] for i in 1:2]
    end

    function fit_modes(observations; u0 = base_u0)
        datasets = [
            (
                x = times,
                y = observations[i],
                state_index = i,
                residual_scale = maximum(observations[i]),
            ) for i in 1:2
        ]
        shared! = function (du, u, p, t)
            for i in eachindex(u)
                du[i] = p[1] * u[i]
            end
        end
        partial! = function (du, u, p, t)
            du[1] = p[1] * exp(-p[2]) * u[1]
            du[2] = p[1] * exp(p[2]) * u[2]
        end
        independent! = function (du, u, p, t)
            du[1] = p[1] * u[1]
            du[2] = p[2] * u[2]
        end
        shared = run_joint_fit(shared!, datasets, u0, [0.2]; bounds = [(0.01, 1.0)], optimizer = :nelder_mead, maxiters = 800)
        partial = run_joint_fit(partial!, datasets, u0, [0.2, 0.0]; bounds = [(0.01, 1.0), (-log(1.05), log(1.05))], optimizer = :nelder_mead, maxiters = 1200)
        independent = run_joint_fit(independent!, datasets, u0, [0.2, 0.2]; bounds = [(0.01, 1.0), (0.01, 1.0)], optimizer = :nelder_mead, maxiters = 1200)
        return (shared = shared, partial = partial, independent = independent)
    end

    identical = fit_modes(generate_pair(0.25, 0.0))
    @test identical.shared.bic < identical.partial.bic

    truth_contrast = log(1.035)
    modest = fit_modes(generate_pair(0.25, truth_contrast))
    @test modest.partial.bic < modest.shared.bic
    @test abs(modest.partial.params[2] - truth_contrast) < 0.01
    @test abs(modest.partial.params[2]) <= log(1.05) + 1e-8

    excessive = fit_modes(generate_pair(0.25, log(1.12)))
    excessive_summary = summarize_pooling_bic(DataFrame(
        cell_line = fill("synthetic", 3),
        model = fill("exponential", 3),
        pooling_mode = ["shared", "partial_5pct", "independent_diagnostic"],
        bic = [excessive.shared.bic, excessive.partial.bic, excessive.independent.bic],
        eligible_for_inheritance = [true, true, false],
    ))
    @test excessive_summary.status.inadequate_pooling[1]

    reference_observations = generate_pair(0.25, 0.0)
    reference = fit_modes(reference_observations).shared
    scaled_observations = [reference_observations[1], 10 .* reference_observations[2]]
    scaled = fit_modes(scaled_observations; u0 = [base_u0[1], 10 * base_u0[2]]).shared
    @test abs(reference.params[1] - scaled.params[1]) < 1e-3
end
