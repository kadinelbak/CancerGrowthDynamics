using CSV
using DataFrames

@testset "Staged report model and identifiability audit" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    report_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(start = package_root)
    html = read(report_path, String)
    teaching_html = first(split(html, "id=\"model-appendix\""))

    @test count("class=\"bic-figure\"", html) == 7
    @test !occursin("<th>BIC</th>", teaching_html)
    @test occursin("<th>delta_BIC</th>", teaching_html)
    @test count("Simplest candidate", teaching_html) >= 4
    @test occursin("Model shown: asymmetric competition with lineage-specific loss", html)
    @test occursin("Model shown: Load scaling plus tolerant-state shift", html)
    @test count("same y-axis range", html) == 5

    @test MechanicalAutomaticModeling.FitWorkflows._shared_y_limits(
        DataFrame(observed = [0.0, 10.0], predicted = [2.0, 20.0]),
    ) ≈ (0.0, 21.0)

    @test occursin("Complete model, BIC, and parameter audit", html)
    appendix_html = last(split(html, "id=\"model-appendix\""))
    @test occursin("<th>BIC</th>", appendix_html)
    @test occursin("<th>Parameters</th>", appendix_html)
    @test occursin("Strict drug-effect inheritance", appendix_html)

    @test occursin("GrowthParameterEstimation two-sided sensitivity and bound profile", html)
    @test occursin("endpoint bootstrap: 95% confidence intervals", html)
    endpoint_path = joinpath(package_root, "outputs", "csv", "coculture_treated", "linked_treatment_endpoint_bootstrap.csv")
    @test isfile(endpoint_path)
    endpoint = CSV.read(endpoint_path, DataFrame)
    @test nrow(endpoint) == 12
    @test all(endpoint.n_wells .== 6)
    @test all(endpoint.ci95_lower .<= endpoint.observed_mean .<= endpoint.ci95_upper)
    @test all(isfinite, endpoint.model_prediction)
    @test all(endpoint.n_bootstrap .== 5000)
end
