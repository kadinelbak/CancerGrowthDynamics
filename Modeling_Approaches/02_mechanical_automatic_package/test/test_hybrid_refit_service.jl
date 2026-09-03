using Test

@testset "Hybrid conditional-refit contract" begin
    service = MechanicalAutomaticModeling.HybridRefitService
    package_root = normpath(joinpath(@__DIR__, ".."))
    request = (
        selections = (
            stage1_A2780Naive = (model = "baranyi_theta_logistic_growth", pooling = "partial_5pct"),
            stage1_A2780cis = (model = "adaptation_theta_logistic_growth", pooling = "shared"),
            stage2_A2780Naive = (model = "joint_ic_effect_hill_ramp_onset", pooling = "partial_5pct"),
            stage2_A2780cis = (model = "joint_ic_effect_two_population", pooling = "shared"),
            stage3 = (model = "lv_asymmetric_competition_death", pooling = "partial_5pct"),
            stage4 = (model = "load_plus_tolerant_growth_context", pooling = "linked_global"),
        ),
        max_time_per_fit = 12.0,
    )
    canonical = service.canonical_refit_request(request)
    @test canonical.schema_version == 1
    @test service.refit_request_id(request) == service.refit_request_id(canonical)
    @test length(service.refit_request_id(request)) == 20
    @test service.validate_refit_request(request; start = package_root).stage4_compatible

    incompatible = merge(request, (selections = merge(request.selections, (
        stage2_A2780cis = (model = "joint_ic_effect_transit_death", pooling = "shared"),
    )),))
    @test !service.validate_refit_request(incompatible; start = package_root).stage4_compatible

    diagnostic = merge(request, (selections = merge(request.selections, (
        stage1_A2780cis = (model = "logistic_growth", pooling = "independent_diagnostic"),
    )),))
    @test_throws ErrorException service.validate_refit_request(diagnostic; start = package_root)
end
