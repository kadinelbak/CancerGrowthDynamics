@testset "Adaptive simulator report contract" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    report_path = joinpath(package_root, "outputs", "reports", "a2780_adaptive_treatment_simulations.html")
    config_path = joinpath(package_root, "outputs", "reports", "a2780_adaptive_simulator_config.js")
    report = read(report_path, String)

    @test isfile(config_path)
    @test occursin("a2780_adaptive_simulator_config.js", report)
    @test occursin("id=\"advancedModelSelect\"", report)
    @test occursin("id=\"showLatent\"", report)
    @test occursin("id=\"addInterval\"", report)
    @test occursin("Active windows cannot overlap", report)
    @test occursin("80% scenario ensemble after day 14", report)
    @test occursin("table-scroll", report)
    @test occursin("doseToEffect", report)
    @test occursin("selectedAdvancedModel = adaptiveConfig.winners.stage4.model", report)
end
