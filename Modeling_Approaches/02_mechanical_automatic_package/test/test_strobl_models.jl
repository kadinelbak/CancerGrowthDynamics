@testset "Strobl literature benchmark models" begin
    workflows = MechanicalAutomaticModeling.FitWorkflows
    variants = workflows.STROBL_MODEL_VARIANTS
    @test length(variants) == 7
    @test length(unique(variant.name for variant in variants)) == 7
    @test Set(variant.paper_model for variant in variants) == Set([
        "2params_noCost_noTurnover",
        "3params_noTurnover",
        "3params_noCost",
        "4params",
        "supplement_cost_in_KR",
        "supplement_cost_in_dR",
        "supplement_density_dependent_death",
    ])

    baseline = (r = 0.2,)
    no_cost = workflows._strobl_variant("strobl_birth_no_cost_no_turnover")
    values = workflows._strobl_parameter_values([100.0, 120.0], [:K_20k, :K_30k], no_cost, "20k", baseline)
    dS, dR = workflows._strobl_rhs(10.0, 10.0, values, 0.0, false)
    @test dS ≈ 1.6
    @test dR ≈ 1.6

    growth_cost = workflows._strobl_variant("strobl_birth_growth_cost_turnover")
    cost_names = [:K_20k, :K_30k, :resistance_cost, :turnover_fraction, :d_D]
    cost_values = workflows._strobl_parameter_values([100.0, 120.0, 0.25, 0.1, 0.75], cost_names, growth_cost, "20k", baseline)
    treated_dS, treated_dR = workflows._strobl_rhs(10.0, 10.0, cost_values, 0.5, false)
    untreated_dS, untreated_dR = workflows._strobl_rhs(10.0, 10.0, cost_values, 0.0, false)
    @test treated_dS < untreated_dS
    @test treated_dR == untreated_dR
    @test cost_values.r_resistant ≈ 0.75 * cost_values.r_sensitive
    @test cost_values.d_sensitive ≈ cost_values.d_resistant

    reduced_names = cost_names[1:4]
    capacity_cost = workflows._strobl_variant("strobl_birth_capacity_cost_turnover")
    capacity_values = workflows._strobl_parameter_values([100.0, 120.0, 0.25, 0.1], reduced_names, capacity_cost, "20k", baseline)
    @test capacity_values.K_resistant ≈ 75.0
    @test capacity_values.r_resistant ≈ capacity_values.r_sensitive

    death_cost = workflows._strobl_variant("strobl_birth_death_cost_turnover")
    death_values = workflows._strobl_parameter_values([100.0, 120.0, 0.25, 0.1], reduced_names, death_cost, "20k", baseline)
    @test death_values.d_resistant ≈ 1.25 * death_values.d_sensitive

    density_death = workflows._strobl_variant("strobl_density_dependent_death")
    density_values = workflows._strobl_parameter_values([100.0, 120.0, 0.25, 0.1], reduced_names, density_death, "20k", baseline)
    density_dS, density_dR = workflows._strobl_rhs(10.0, 10.0, density_values, 0.0, true)
    @test isfinite(density_dS)
    @test isfinite(density_dR)
end
