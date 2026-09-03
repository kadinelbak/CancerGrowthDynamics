using Test

@testset "Conditional model-path tournament" begin
    tournament = MechanicalAutomaticModeling.ModelPathTournament
    package_root = normpath(joinpath(@__DIR__, ".."))
    config = tournament.TournamentConfig(
        candidates_per_lineage = 2,
        beam_width = 2,
        multistarts = 1,
        max_time_per_fit = 0.1,
    )
    beam = tournament.stage1_beam(; start = package_root, config = config)
    @test size(beam, 1) >= 2
    @test issorted(beam.stage1_bic)
    @test all(isfinite, beam.stage1_bic)
    @test any(beam.simple_path)

    section = tournament.tournament_report_section(joinpath(package_root, "outputs", "model_path_tournament"))
    @test contains(section, "Conditional Model Tournament")
    @test contains(section, "Stage 2 A2780Naive")
    @test contains(section, "Blocked validation")
end
