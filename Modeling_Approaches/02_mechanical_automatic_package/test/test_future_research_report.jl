@testset "Future research report contract" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    report_path = joinpath(package_root, "outputs", "reports", "a2780_future_research.html")
    docs_report_path = normpath(joinpath(package_root, "..", "..", "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports", "a2780_future_research.html"))

    @test isfile(report_path)
    @test isfile(docs_report_path)

    report = read(report_path, String)
    docs_report = read(docs_report_path, String)
    @test report == docs_report
    @test occursin("A2780 Coculture Adaptive Therapy", report)
    @test occursin("<summary>Contents</summary>", report)
    @test occursin("id=\"core-model\"", report)
    @test !occursin("The Smallest Defensible Model", report)
    @test !occursin("model-flow", report)
    @test occursin("D_{cmd}(t)", report)
    @test occursin("D_{k+1}", report)
    @test occursin("Policy choices", report)
    @test occursin("When adaptive therapy can work", report)
    @test occursin("Low resources: start with one fitted effect", report)
    @test occursin("Low resources: separate lineage effects only if needed", report)
    @test occursin("Cisplatin: start with direct on/off treatment", report)
    @test occursin("Cisplatin: add delayed exposure only if needed", report)
    @test occursin("\\beta_N\\,u(t)\\,N", report)
    @test occursin("\\frac{dE}{dt}=k_E[u(t)-E]", report)
    @test occursin("the one fitted parameter", report)
    @test occursin("\\rho_N", report)
    @test !occursin("q_N(E)", report)
    @test !occursin("G_N(N", report)
    @test !occursin("G_C(C", report)
    @test !occursin("L_N", report)
    @test !occursin("L_C", report)
    @test occursin("id=\"compactRanking\"", report)
    @test occursin("id=\"paperRanking\"", report)
    @test occursin("class=\"annotated-bibliography\"", report)
    @test occursin("<summary>Annotated bibliography <span class=\"small\">15 papers</span></summary>", report)
    @test !occursin("id=\"roadmap\"", report)
    @test !occursin("Future Research Roadmap", report)
    @test occursin("function buildCompactRanking()", report)
    @test occursin("Biological relevance", report)
    @test occursin("Annotated bibliography", report)
    @test !occursin("The 15 Most Pertinent Papers", report)
    @test !occursin("implementation horizons", report)
    @test occursin("The Minimum Experimental Program", report)
    @test occursin("Interpretive boundary", report)
    @test occursin("Not peer reviewed", report)
    @test occursin("A larger carrying capacity alone cannot represent patient tumors", report)
    @test occursin("keep latent cis states out of the default graph", report)

    ranked_rows = collect(eachmatch(r"<td class=\"rank\">(\d+)</td>", report))
    @test length(ranked_rows) == 15
    @test parse.(Int, getindex.(getfield.(ranked_rows, :captures), 1)) == collect(1:15)

    for doi in (
        "10.1158/0008-5472.CAN-25-0351",
        "10.1016/S1476-5586(03)80008-8",
        "10.1016/j.cels.2024.04.003",
        "10.3389/fonc.2024.1304691",
        "10.1371/journal.pcbi.1012073",
        "10.1158/0008-5472.CAN-20-0806",
        "10.1158/0008-5472.CAN-08-3658",
    )
        @test occursin(doi, report)
    end
end
