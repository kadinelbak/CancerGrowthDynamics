using Test
using DifferentialEquations
using GrowthParameterEstimation
using DataFrames
using CSV
using Random

@testset "GrowthParameterEstimation" begin
    @test isdefined(GrowthParameterEstimation, :Models)

    # Basic solve of logistic growth
    u0 = [1.0]
    p  = (0.5, 10.0)  # tuple avoids param-scalar inference
    tspan = (0.0, 2.0)
    prob = ODEProblem(GrowthParameterEstimation.Models.logistic_growth!, u0, tspan, p)
    sol  = solve(prob, Tsit5(); saveat = 0.5)
    @test sol.u[end][1] > 0  # check scalar state value, not the whole vector

    # Smoke test: evaluate BIC on a known solution
    x = [0.0, 1.0, 2.0, 3.0]
    y = [1.0, 1.8, 2.6, 3.4]
    prob_fit = ODEProblem(GrowthParameterEstimation.Models.logistic_growth!, [y[1]], (x[1], x[end]), [0.1, 5.0])
    bic, ssr = GrowthParameterEstimation.Fitting.calculate_bic(prob_fit, x, y, Tsit5(), [0.1, 5.0])
    @test ssr >= 0
    @test isfinite(bic)

    # Smoke test fitting routine (ensures run_single_fit works end-to-end)
    fit = GrowthParameterEstimation.run_single_fit(x, y, [0.1, 5.0]; solver = Tsit5(), max_time = 1.0, show_stats = false)
    @test length(fit.params) == 2
    @test isfinite(fit.bic)

    @testset "Data Schema Validation" begin
        raw = DataFrame(time = [0.0, 1.0, 2.0], count = [1.0, 2.0, 3.1], dose = [0.0, 0.0, 0.0], cell_line = ["A", "A", "A"], density = [1.0, 1.0, 1.0], replicate = [1, 1, 1])
        norm = normalize_schema(raw)
        @test validate_timeseries(norm)
        @test all(c -> c in Symbol.(names(norm)), REQUIRED_COLUMNS)
        @test all(norm.treatment_amount .== norm.dose)

        raw_treatment = DataFrame(time = [0.0, 1.0], count = [2.0, 3.0], treatment_amount = [1.5, 1.5], cell_line = ["A", "A"], density = [1.0, 1.0], replicate = [1, 1])
        norm_treatment = normalize_schema(raw_treatment)
        @test all(norm_treatment.dose .== 1.5)
        @test all(norm_treatment.treatment_amount .== 1.5)

        csv_path = joinpath(tempdir(), "gpe_test_timeseries.csv")
        CSV.write(csv_path, norm)
        loaded = load_timeseries(csv_path)
        @test nrow(loaded) == nrow(norm)
    end

    @testset "Exposure + Observation" begin
        times = collect(0.0:1.0:5.0)
        exp_constant = build_exposure(:constant; value = 0.25)
        exp_pulse = build_exposure(:pulse; amplitude = 1.0, start_time = 2.0, end_time = 4.0)
        exp_step = build_exposure(:stepped; change_times = [0.0, 3.0], values = [0.1, 0.4])
        exp_decay = build_exposure(:decay; c0 = 1.0, decay_rate = 0.3, t0 = 1.0)

        @test length(evaluate_exposure(exp_constant, times)) == length(times)
        @test evaluate_exposure(exp_pulse, [1.0, 3.0]) == [0.0, 1.0]
        @test length(evaluate_exposure(exp_step, times)) == length(times)
        @test length(evaluate_exposure(exp_decay, times)) == length(times)

        obs = ObservationSpec("viable_scaled", viable_total, 1.2, 0.1)
        state = [3.0, 2.0]
        @test isapprox(observed_signal(obs, state, nothing, 0.0), 3.6)
        sum_fn = sum_states([1, 2])
        @test sum_fn(state, nothing, 0.0) == 5.0
    end

    @testset "Model RHS Smoke" begin
        du = [0.0]
        logistic_growth!(du, [2.0], [0.2, 10.0], 0.0)
        @test isfinite(du[1])

        logistic_growth_with_death!(du, [2.0], [0.2, 10.0, 0.01], 0.0)
        @test isfinite(du[1])

        gompertz_growth!(du, [2.0], [0.2, 1.0, 10.0], 0.0)
        @test isfinite(du[1])

        gompertz_growth_with_death!(du, [2.0], [0.2, 1.0, 10.0, 0.01], 0.0)
        @test isfinite(du[1])

        exponential_growth!(du, [2.0], [0.2], 0.0)
        @test isfinite(du[1])

        exponential_growth_with_delay!(du, [2.0], [0.2, 10.0, 1.0], 2.0)
        @test isfinite(du[1])

        logistic_growth_with_delay!(du, [2.0], [0.2, 10.0, 1.0], 2.0)
        @test isfinite(du[1])

        exponential_growth_with_death_and_delay!(du, [2.0], [0.2, 10.0, 0.01, 1.0], 2.0)
        @test isfinite(du[1])
    end

    @testset "Registry + Simulation" begin
        models = list_models()
        @test "logistic_growth" in models
        @test "null_coculture" in models
        @test "lotka_volterra_competition" in models
        @test "lotka_volterra_hill_competition" in models

        spec = get_model("theta_logistic_hill_inhibition")
        times = collect(0.0:1.0:5.0)
        params = [0.4, 50.0, 1.0, 0.5, 1.5]
        sim = simulate(spec, times, params; u0 = [1.0], exposure = ConstantExposure(0.2))
        @test sim.success
        @test length(sim.observed) == length(times)
        @test all(isfinite, sim.observed)

        sweep = run_sweep(
            get_model("sensitive_resistant"),
            [0.4, 0.2, 100.0, 0.05, 0.7, 0.5, 2.0],
            SweepGrid([10.0, 20.0], [0.1, 0.3], [0.0, 1.0], collect(0.0:1.0:4.0)),
        )
        @test nrow(sweep.summary) == 8
        @test all(col -> col in Symbol.(names(sweep.summary)), [:seed_total, :resistant_fraction, :dose, :final_total])
    end

    @testset "Pipeline Fit + Ranking" begin
        df = DataFrame(
            time = vcat(collect(0.0:1.0:6.0), collect(0.0:1.0:6.0)),
            count = vcat([1.0, 1.5, 2.2, 3.0, 3.8, 4.5, 5.1], [1.0, 1.3, 1.9, 2.5, 3.1, 3.6, 4.0]),
            error = fill(0.2, 14),
            dose = vcat(fill(0.2, 7), fill(0.6, 7)),
            cell_line = fill("A549", 14),
            density = fill(1.0, 14),
            replicate = vcat(fill(1, 7), fill(2, 7)),
            unit_time = fill("h", 14),
            unit_count = fill("count", 14),
        )

        conditions = build_conditions(df)
        ranked = rank_models(["logistic_growth", "theta_logistic_hill_inhibition"], conditions; n_starts = 3, maxiters = 80, top_k = 2, seed = 7)

        @test nrow(ranked.ranking) == 2
        @test ranked.ranking.bic[1] <= ranked.ranking.bic[2]
        finite_deltas = filter(isfinite, ranked.ranking.delta_bic)
        @test all(finite_deltas .>= 0)

        # Workflow config + plotting + export + pipeline run
        cfg = default_config(output_dir = joinpath(tempdir(), "gpe_pipeline_test"))
        cfg_path = save_config(joinpath(tempdir(), "gpe_cfg_test.toml"), cfg)
        cfg_loaded = load_config(cfg_path)
        @test cfg_loaded.top_k == cfg.top_k

        # Use a single-condition subset for deterministic export behavior
        df_one = df[df.replicate .== 1, :]
        conditions_one = build_conditions(df_one)
        ranked_logistic = rank_models(["logistic_growth"], conditions_one; n_starts = 3, maxiters = 80, top_k = 1, seed = 11)
        generated = plot_topk(ranked_logistic; conditions = conditions_one, top_k = 1, output_dir = joinpath(tempdir(), "gpe_plot_topk_test"))
        @test generated isa Vector{String}

        try
            exports = export_results(ranked_logistic; output_dir = joinpath(tempdir(), "gpe_export_test"))
            @test isfile(exports.ranking)
            @test isfile(exports.params)
            @test isfile(exports.failures)
        catch err
            @test occursin("No successful model fits available to export", string(err))
        end

        pipe = run_pipeline(df; config = cfg_loaded, include_models = ["logistic_growth"])
        @test nrow(pipe.ranking) == 1

        df_stage = DataFrame(
            time = vcat(collect(0.0:1.0:5.0), collect(0.0:1.0:5.0)),
            count = vcat([1.0, 1.6, 2.4, 3.1, 3.7, 4.2], [1.0, 1.3, 1.7, 2.0, 2.3, 2.5]),
            error = fill(0.15, 12),
            dose = vcat(fill(0.0, 6), fill(0.8, 6)),
            cell_line = fill("A549", 12),
            density = fill(1.0, 12),
            replicate = vcat(fill(1, 6), fill(1, 6)),
            unit_time = fill("h", 12),
            unit_count = fill("count", 12),
            culture_type = fill("monoculture", 12),
            population_type = fill("naive", 12),
        )

        stages = [
            PipelineStage(
                "untreated_monoculture",
                "Untreated mono",
                row -> row[:culture_type] == "monoculture" && row[:dose] == 0.0,
                [:cell_line, :population_type, :replicate],
                ["logistic_growth"],
                Symbol[],
                Dict{Symbol,Float64}(),
                Dict{Symbol,Tuple{String,Symbol}}(),
            ),
            PipelineStage(
                "treated_monoculture",
                "Treated mono",
                row -> row[:culture_type] == "monoculture" && row[:dose] > 0.0,
                [:dose, :cell_line, :population_type, :replicate],
                ["theta_logistic_hill_inhibition"],
                [:ic50, :hill],
                Dict{Symbol,Float64}(:theta => 1.0),
                Dict(
                    :r => ("untreated_monoculture", :r),
                    :K => ("untreated_monoculture", :K),
                ),
            ),
        ]

        staged = run_staged_pipeline(df_stage; stages=stages, config=cfg_loaded, selection_mode=:best_bic, export_stage_results=false)
        @test staged.completed
        @test length(staged.stages) == 2
        @test staged.stages[1].selected_model == "logistic_growth"
        @test haskey(staged.parameter_bank, "untreated_monoculture")

        staged_checkpoint = run_staged_pipeline(df_stage; stages=stages, config=cfg_loaded, selection_mode=:manual, manual_choices=Dict("untreated_monoculture" => "logistic_growth"), export_stage_results=false)
        @test staged_checkpoint.halted_stage == "treated_monoculture"
        @test staged_checkpoint.stages[end].status == "awaiting_selection"

        pop_stages = default_population_stages(["naive", "resistant"])
        @test any(s -> s.name == "untreated_monoculture_naive", pop_stages)
        @test any(s -> s.name == "treated_monoculture_resistant", pop_stages)

        df_cellline = DataFrame(
            time = vcat(collect(0.0:1.0:3.0), collect(0.0:1.0:3.0)),
            count = vcat([1.0, 1.3, 1.8, 2.2], [0.9, 1.2, 1.5, 1.8]),
            error = fill(0.1, 8),
            treatment_amount = vcat(fill(0.0, 4), fill(1.0, 4)),
            cell_line = vcat(fill("A549", 4), fill("H1975", 4)),
            density = fill(1.0, 8),
            replicate = fill(1, 8),
            culture_type = fill("monoculture", 8),
            population_type = fill("naive", 8),
        )

        cl_stages = default_population_cellline_stages(df_cellline; populations=["naive", "resistant"])
        @test any(s -> s.name == "treated_monoculture_naive_a549", cl_stages)
        @test any(s -> s.name == "treated_monoculture_naive_h1975", cl_stages)
        @test any(s -> s.name == "treated_coculture_a549", cl_stages)

        st_a549 = only(filter(s -> s.name == "treated_coculture_a549", cl_stages))
        @test st_a549.inherited_params[:ic50S] == ("treated_monoculture_naive_a549", :ic50)
        @test st_a549.inherited_params[:ic50R] == ("treated_monoculture_resistant_a549", :ic50)
        @test st_a549.inherited_params[:rS] == ("untreated_monoculture_naive", :r)
        @test st_a549.inherited_params[:rR] == ("untreated_monoculture_resistant", :r)

        df_summary = copy(df_stage)
        df_summary.culture_type .= "monoculture"
        df_summary.population_type .= "naive"
        df_summary.treatment_amount = df_summary.dose
        df_summary.ic50_reference = fill(2.2, nrow(df_summary))
        summary_table = summarize_datasets(df_summary)
        @test nrow(summary_table) >= 1
        @test all(col -> col in Symbol.(names(summary_table)), [:duration_days, :treatment_amount, :ic50_reference, :n_timepoints])
        @test all(summary_table.duration_days .>= 0)

        strict_df = copy(df_summary)
        @test validate_strict_schema(strict_df)

        strict_bad = select(strict_df, Not(:culture_type))
        @test_throws ErrorException validate_strict_schema(strict_bad)

        qc = generate_qc_report(strict_df)
        @test all(col -> col in Symbol.(names(qc.missingness)), [:column, :n_missing, :frac_missing])
        qc_paths = save_qc_report(qc; output_dir=joinpath(tempdir(), "gpe_qc_test"))
        @test isfile(qc_paths.missingness)
        @test isfile(qc_paths.condition_summary)

        boot_spec = get_model("logistic_growth")
        boot_unc = bootstrap_stage_uncertainty(
            boot_spec,
            df_stage[df_stage.dose .== 0.0, :];
            condition_cols=[:cell_line, :population_type, :replicate],
            shared_params=[:r, :K],
            n_bootstrap=3,
            n_starts=2,
            maxiters=60,
            seed=21,
        )
        @test haskey(boot_unc, :r)
        @test haskey(boot_unc[:r], "ci_lower")

        cfg_resume = default_config(output_dir = joinpath(tempdir(), "gpe_resume_test"))
        staged_resume_source = run_staged_pipeline(
            df_stage;
            stages=stages,
            config=cfg_resume,
            selection_mode=:manual,
            manual_choices=Dict("untreated_monoculture" => "logistic_growth"),
            export_stage_results=true,
            strict_schema=true,
            n_bootstrap=2,
        )
        @test isfile(staged_resume_source.manifest_path)

        manifest = load_run_manifest(staged_resume_source.manifest_path)
        @test haskey(manifest.parameter_bank, "untreated_monoculture")

        staged_resumed = run_staged_pipeline(
            df_stage;
            stages=stages,
            config=cfg_resume,
            selection_mode=:manual,
            manual_choices=Dict(
                "untreated_monoculture" => "logistic_growth",
                "treated_monoculture" => "theta_logistic_hill_inhibition",
            ),
            export_stage_results=false,
            resume_manifest_path=staged_resume_source.manifest_path,
            resume_from_stage="treated_monoculture",
        )
        @test staged_resumed.completed
        @test haskey(staged_resumed.parameter_bank, "untreated_monoculture")
    end

    @testset "Fitting Compare APIs" begin
        x = collect(0.0:1.0:6.0)
        y = [1.0, 1.5, 2.1, 2.9, 3.7, 4.3, 4.9]

        comp = compare_models(
            x, y,
            "Logistic", logistic_growth!, [0.2, 8.0],
            "Gompertz", gompertz_growth!, [0.2, 1.0, 8.0];
            solver = Tsit5(),
            bounds1 = [(0.01, 2.0), (2.0, 100.0)],
            bounds2 = [(0.01, 2.0), (0.1, 5.0), (2.0, 100.0)],
            show_stats = false,
            output_csv = joinpath(tempdir(), "gpe_compare_models_test.csv"),
        )
        @test haskey(comp, :best_model)

        compare_datasets(
            x, y, "A", logistic_growth!, [0.2, 8.0],
            x, y .* 0.95, "B", logistic_growth!, [0.2, 8.0];
            solver = Tsit5(),
            bounds1 = [(0.01, 2.0), (2.0, 100.0)],
            bounds2 = [(0.01, 2.0), (2.0, 100.0)],
            show_stats = false,
            output_csv = joinpath(tempdir(), "gpe_compare_datasets_test.csv"),
        )
        @test isfile(joinpath(tempdir(), "gpe_compare_datasets_test.csv"))

        specs = Dict(
            "Logistic" => (model = logistic_growth!, p0 = [0.2, 8.0], bounds = [(0.01, 2.0), (2.0, 100.0)]),
            "Exp" => (model = exponential_growth!, p0 = [0.2], bounds = [(0.001, 2.0)]),
        )
        dict_fits = compare_models_dict(x, y, specs; default_solver = Tsit5(), show_stats = false, output_csv = joinpath(tempdir(), "gpe_compare_dict_test.csv"))
        @test haskey(dict_fits, "Logistic")

        three_named = fit_three_datasets(
            x, y, "A",
            x, y .* 0.9, "B",
            x, y .* 1.1, "C",
            [0.2, 8.0];
            model = logistic_growth!,
            solver = Tsit5(),
            bounds = [(0.01, 2.0), (2.0, 100.0)],
            show_stats = false,
            output_csv = joinpath(tempdir(), "gpe_three_named_test.csv"),
        )
        @test haskey(three_named, :fit1)

        x_many = [collect(x), collect(x)]
        y_many = [collect(y), collect(y .* 0.95)]
        three_many = fit_three_datasets(x_many, y_many; model = logistic_growth!, solver = Tsit5(), bounds = [(0.01, 2.0), (2.0, 100.0)])
        @test three_many.summary.n_total == 2
    end

    @testset "Analysis APIs" begin
        Random.seed!(7)
        x = collect(0.0:1.0:8.0)
        y = [1.0, 1.6, 2.4, 3.1, 3.9, 4.6, 5.2, 5.7, 6.0]
        p0 = [0.2, 12.0]
        bounds = [(0.01, 2.0), (2.0, 100.0)]

        fit = run_single_fit(x, y, p0; model = logistic_growth!, solver = Tsit5(), bounds = bounds, show_stats = false)

        loo = leave_one_out_validation(x, y, p0; model = logistic_growth!, solver = Tsit5(), bounds = bounds, show_stats = false)
        @test isfinite(loo.rmse)

        kfold = k_fold_cross_validation(x, y, p0; k_folds = 3, model = logistic_growth!, solver = Tsit5(), bounds = bounds, show_stats = false)
        @test isfinite(kfold.overall_rmse)

        sens = parameter_sensitivity_analysis(x, y, fit; perturbation = 0.1, model = logistic_growth!, solver = Tsit5())
        @test length(sens.ranking) >= 1

        resid = residual_analysis(x, y, fit; model = logistic_growth!, solver = Tsit5())
        @test isfinite(resid.statistics.rmse)

        enh = enhanced_bic_analysis(
            x,
            y;
            models = [logistic_growth!],
            model_names = ["Logistic"],
            p0_values = [[0.2, 12.0]],
            solver = Tsit5(),
        )
        @test enh.best_model.model_name == "Logistic"
    end

    @testset "Joint Fitting APIs" begin
        function logistic_joint!(du, u, p, t)
            r, K = p
            du[1] = r * u[1] * (1 - u[1] / K)
            du[2] = r * u[2] * (1 - u[2] / K)
        end

        x1 = collect(0.0:1.0:5.0)
        x2 = collect(0.0:1.0:5.0)
        u0_joint = [1.0, 2.0]
        p_true = [0.35, 25.0]
        prob_joint = ODEProblem(logistic_joint!, u0_joint, (x1[1], x1[end]), p_true)
        sol_joint = solve(prob_joint, Tsit5(); saveat = x1)
        y1 = [u[1] for u in sol_joint.u]
        y2 = [u[2] for u in sol_joint.u]

        datasets = [
            (x = x1, y = y1, state_index = 1),
            (x = x2, y = y2, state_index = 2),
        ]

        joint_fit = run_joint_fit(logistic_joint!, datasets, u0_joint, [0.2, 20.0]; bounds = [(0.01, 1.0), (5.0, 60.0)], show_stats = false)
        @test length(joint_fit.params) == 2
        @test isfinite(joint_fit.bic)
        @test joint_fit.sse >= 0

        delayed_times = collect(1.0:1.0:5.0)
        delayed_prob = ODEProblem(logistic_joint!, u0_joint, (0.0, delayed_times[end]), p_true)
        delayed_sol = solve(delayed_prob, Tsit5(); saveat = delayed_times)
        delayed_data = [
            (x = delayed_times, y = [u[1] for u in delayed_sol.u], state_index = 1),
            (x = delayed_times, y = [u[2] for u in delayed_sol.u], state_index = 2),
        ]
        delayed_fit = run_joint_fit(
            logistic_joint!, delayed_data, u0_joint, [0.2, 20.0];
            bounds = [(0.01, 1.0), (5.0, 60.0)],
            initial_time = 0.0,
            show_stats = false,
        )
        @test delayed_fit.initial_time == 0.0
        @test first(delayed_fit.save_times) == 1.0
        @test all(length(prediction) == length(delayed_times) for prediction in delayed_fit.predictions)
        @test_throws ArgumentError run_joint_fit(
            logistic_joint!, delayed_data, u0_joint, [0.2, 20.0];
            bounds = [(0.01, 1.0), (5.0, 60.0)],
            initial_time = 2.0,
            show_stats = false,
        )

        scaled_datasets = [merge(ds, (residual_scale = maximum(ds.y),)) for ds in datasets]
        scaled_joint_fit = run_joint_fit(
            logistic_joint!, scaled_datasets, u0_joint, [0.2, 20.0];
            bounds = [(0.01, 1.0), (5.0, 60.0)],
            show_stats = false,
        )
        @test isfinite(scaled_joint_fit.bic)
        @test scaled_joint_fit.scaled_sse == scaled_joint_fit.sse
        @test isfinite(scaled_joint_fit.raw_sse)
        @test scaled_joint_fit.raw_sse >= scaled_joint_fit.scaled_sse

        function split_joint!(du, u, p, t)
            decay_live, decay_damaged, _ = p
            du[1] = -decay_live * u[1]
            du[2] = decay_live * u[1] - decay_damaged * u[2]
        end
        split_times = collect(0.0:1.0:4.0)
        split_true = [0.4, 0.2, 0.75]
        split_u0_builder(p) = [10.0 * p[3], 10.0 * (1 - p[3])]
        split_prob = ODEProblem(split_joint!, split_u0_builder(split_true), (0.0, 4.0), split_true)
        split_sol = solve(split_prob, Tsit5(); saveat = split_times)
        split_y = [u[1] + 0.5 * u[2] for u in split_sol.u]
        split_data = [(
            x = split_times,
            y = split_y,
            observable = (u, p, t) -> u[1] + 0.5 * u[2],
        )]
        split_fit = run_joint_fit(
            split_joint!, split_data, [7.0, 3.0], [0.3, 0.3, 0.6];
            bounds = [(0.01, 1.0), (0.01, 1.0), (0.01, 0.99)],
            u0_builder = split_u0_builder,
            show_stats = false,
        )
        @test length(split_fit.params) == 3
        @test isfinite(split_fit.bic)
        @test all(isfinite, only(split_fit.predictions))

        specs = Dict(
            "LogisticJoint" => (model = logistic_joint!, p0 = [0.2, 20.0], bounds = [(0.01, 1.0), (5.0, 60.0)]),
        )
        joint_models = compare_joint_models_dict(datasets, u0_joint, specs; default_solver = Tsit5(), show_stats = false, output_csv = joinpath(tempdir(), "joint_compare_smoke.csv"))
        @test haskey(joint_models, "LogisticJoint")

        multistart = run_joint_multistart(
            logistic_joint!,
            scaled_datasets,
            u0_joint,
            [[0.08, 12.0], [0.3, 30.0]];
            bounds = [(0.01, 1.0), (5.0, 60.0)],
            optimizer = :nelder_mead,
            maxiters = 800,
            show_stats = false,
        )
        @test isfinite(multistart.fit.bic)
        @test nrow(multistart.summary) == 2
        @test multistart.best_start_index in 1:2
        @test all(multistart.summary.status .== "completed")
    end
end

include("joint_summary_test.jl")
include("pooling_profile_test.jl")
