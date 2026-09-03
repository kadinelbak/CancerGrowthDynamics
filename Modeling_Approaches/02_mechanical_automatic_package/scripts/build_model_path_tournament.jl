using MechanicalAutomaticModeling

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
config = MechanicalAutomaticModeling.ModelPathTournament.TournamentConfig(
    candidates_per_lineage = parse(Int, get(ENV, "A2780_TOURNAMENT_STAGE1", "3")),
    beam_width = parse(Int, get(ENV, "A2780_TOURNAMENT_BEAM", "5")),
    simple_delta_bic = parse(Float64, get(ENV, "A2780_TOURNAMENT_SIMPLE_DELTA", "10")),
    multistarts = parse(Int, get(ENV, "A2780_MULTISTARTS", "20")),
    bootstrap_replicates = parse(Int, get(ENV, "A2780_TOURNAMENT_BOOTSTRAPS", "200")),
    blocked_validation = lowercase(get(ENV, "A2780_TOURNAMENT_BLOCKED_VALIDATION", "true")) == "true",
    max_time_per_fit = parse(Float64, get(ENV, "A2780_TOURNAMENT_FIT_SECONDS", "8")),
    force = lowercase(get(ENV, "A2780_TOURNAMENT_FORCE", "false")) == "true",
)
result = MechanicalAutomaticModeling.ModelPathTournament.run_model_path_tournament!(
    start = PACKAGE_ROOT,
    config = config,
)
println("report=", result.report)
println("docs_report=", result.docs_report)
println("ranking=", joinpath(result.output, "complete_path_ranking.csv"))
