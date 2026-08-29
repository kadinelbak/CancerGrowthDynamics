@testset "Reduced-stage comparison report contract" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    report_path = joinpath(package_root, "outputs", "reports", "a2780_reduced_stage_comparison.html")
    docs_report_path = normpath(joinpath(package_root, "..", "..", "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports", "a2780_reduced_stage_comparison.html"))

    @test isfile(report_path)
    @test isfile(docs_report_path)
    report = read(report_path, String)
    @test report == read(docs_report_path, String)
    @test occursin("Reduced-stage treated coculture comparison", report)
    @test occursin("Stage 1 &rarr; Stage 3 &rarr; Stage 4", report)
    @test occursin("Stage 1 &rarr; Stage 2 &rarr; Stage 4", report)
    @test occursin("<em>n</em> = 168", report)
    @test occursin("GrowthParameterEstimation.run_joint_multistart", report)
    @test occursin("Delta BIC", report)
    @test occursin("Every panel uses the same y-axis scale", report)
    @test occursin("Scientific decision", report)

    ranking_path = joinpath(package_root, "outputs", "csv", "reduced_stage_comparison", "reduced_stage_model_ranking.csv")
    ranking = CSV.read(ranking_path, DataFrame)
    @test nrow(ranking) == 17
    @test all(ranking.n_points .== 168)
    @test all(isfinite, ranking.bic)
    @test Set(ranking.route) == Set(["Stage 1 -> Stage 3 -> Stage 4", "Stage 1 -> Stage 2 -> Stage 4"])
    @test minimum(ranking[ranking.route .== "Stage 1 -> Stage 3 -> Stage 4", :bic]) <
          minimum(ranking[ranking.route .== "Stage 1 -> Stage 2 -> Stage 4", :bic])

    image_dir = joinpath(package_root, "outputs", "images", "reduced_stage_comparison")
    @test all(isfile, joinpath.(Ref(image_dir), [
        "competition_first_winner.png",
        "competition_first_simple.png",
        "treatment_first_winner.png",
    ]))
end
