@testset "Untreated lag-growth candidates" begin
    FW = MechanicalAutomaticModeling.FitWorkflows
    initial_counts = [67.0, 100.0]

    hard_lag = FW._untreated_base_spec("lagged_theta_logistic_growth", 4200.0, initial_counts)
    baranyi = FW._untreated_base_spec("baranyi_theta_logistic_growth", 4200.0, initial_counts)
    adaptation = FW._untreated_base_spec("adaptation_theta_logistic_growth", 4200.0, initial_counts)

    @test hard_lag.param_names == [:r, :K, :theta, :lag_time]
    @test baranyi.param_names == [:r, :K, :theta, :q0]
    @test adaptation.param_names == [:r, :K, :theta, :adaptation_rate]

    hard_ode! = FW._untreated_joint_ode("lagged_theta_logistic_growth", hard_lag, "shared")
    du = zeros(2)
    u = copy(initial_counts)
    hard_ode!(du, u, [0.5, 4000.0, 2.0, 2.0], 1.0)
    @test du == [0.0, 0.0]
    hard_ode!(du, u, [0.5, 4000.0, 2.0, 2.0], 3.0)
    @test all(du .> 0)

    baranyi_ode! = FW._untreated_joint_ode("baranyi_theta_logistic_growth", baranyi, "shared")
    baranyi_ode!(du, u, [0.5, 4000.0, 2.0, 0.05], 0.0)
    early_baranyi = copy(du)
    baranyi_ode!(du, u, [0.5, 4000.0, 2.0, 0.05], 8.0)
    @test all(du .> early_baranyi)

    adaptation_ode! = FW._untreated_joint_ode("adaptation_theta_logistic_growth", adaptation, "shared")
    adaptation_ode!(du, u, [0.5, 4000.0, 2.0, 0.7], 0.0)
    @test du == [0.0, 0.0]
    adaptation_ode!(du, u, [0.5, 4000.0, 2.0, 0.7], 3.0)
    @test all(du .> 0)

    pooled = FW._pooling_spec(adaptation, "partial_5pct", ["20k", "30k"])
    @test length(pooled.p0) == 6
    @test pooled.param_names[end-1:end] == [:log_contrast_r, :log_contrast_K]
end
