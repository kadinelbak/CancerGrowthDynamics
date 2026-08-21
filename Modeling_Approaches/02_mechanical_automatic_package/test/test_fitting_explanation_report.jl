@testset "Fitting explanation report" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    repository_root = normpath(joinpath(package_root, "..", ".."))
    canonical = joinpath(package_root, "outputs", "reports", "a2780_fitting_explanation.html")
    mirrored = joinpath(repository_root, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports", "a2780_fitting_explanation.html")

    @test isfile(canonical)
    @test isfile(mirrored)
    @test read(canonical) == read(mirrored)

    html = read(canonical, String)
    @test occursin("Logistic or theta-logistic only", html)
    @test occursin("theta-logistic is the preferred simple baseline for both cell lines", html)
    @test occursin("Downstream stages must be refitted", html)
    @test occursin("a2780_fitting_explanation.html", read(joinpath(repository_root, "docs", "index.html"), String))
end
