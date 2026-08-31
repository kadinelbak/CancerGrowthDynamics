using MechanicalAutomaticModeling

host = get(ENV, "A2780_REFIT_HOST", "127.0.0.1")
port = parse(Int, get(ENV, "A2780_REFIT_PORT", "8766"))
start = get(ENV, "A2780_REPO_START", pwd())

MechanicalAutomaticModeling.HybridRefitService.serve_hybrid_refits(; host, port, start)
