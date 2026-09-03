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
    ) == (0.0, 21.0)

    @test occursin("Complete model, BIC, and parameter audit", html)
    @test occursin("<details class=\"appendix-disclosure\">", html)
    @test occursin("Show appendix", html)
    @test !occursin("<details class=\"appendix-disclosure\" open", html)
    @test occursin("Strobl density-dependent-birth family", html)
    @test occursin("Strobl linked-treatment benchmarks", html)
    @test occursin("Selected for inheritance", html)
    appendix_html = last(split(html, "id=\"model-appendix\""))
    @test occursin("<th>BIC</th>", appendix_html)
    @test occursin("<th>Parameters</th>", appendix_html)
    @test occursin("Strict drug-effect inheritance", appendix_html)

    @test count("Compare fitted assumptions", html) == 4
    @test count("class=\"model-options\"", html) == 6
    @test count("class=\"pooling-options\"", html) == 6
    @test occursin("Choose a model, then its fitted pooling method", html)
    @test occursin("id=\"model-path-data\"", html)
    @test occursin("\"stage2_A2780cis\":\"joint_ic_effect_two_population|shared\"", html)
    @test occursin("logistic_growth|shared", html)
    @test occursin("Selected inheritance path", html)
    @test occursin("Running Delta BIC penalty", html)
    @test occursin("Parameters carried from this fit", html)
    @test occursin("a conditional refit is required", html)
    @test occursin("What pooling means", html)
    @test occursin("Shared pooling", html)
    @test occursin("Partial pooling", html)
    @test occursin("Independent diagnostic", html)
    @test occursin("Linked global", html)
    @test occursin("Well and sample aggregation", html)
    @test occursin("partial_5pct", html)
    @test occursin("Refit downstream", html)
    @test occursin("Conditionally refitted", html)
    @test occursin("http://127.0.0.1:8766", html)
    @test occursin("content-addressed by the exact model and pooling path", html)
    @test occursin("How preview and refitting work", html)
    @test occursin("Partially refitted", html)
    @test count("class=\"live-equation\"", html) == 4
    @test occursin("hot inherited forward preview", html)
    @test occursin("Stage 3 now propagates the selected Stage 1 growth equations and parameters", html)
    @test occursin("effective_parameters", html)
    @test occursin("stage3ForwardRows", html)
    @test !occursin("Its plot remains a stage-local shape comparison", html)

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
