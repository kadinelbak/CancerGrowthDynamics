# Copy the NN-track pages (report.html + figures of every stage, and the index) from outputs/ to the
# published site docs/outputs/nn_track/. CSVs, logs and model files stay in outputs/.
#   julia scripts/publish_docs.jl
const SRC = normpath(joinpath(@__DIR__, "..", "outputs"))
const DST = normpath(joinpath(@__DIR__, "..", "..", "..", "docs", "outputs", "nn_track"))
const STAGES = ["data", "nextday", "ude", "symbolic", "importance", "control", "comparison", "learning_curve", "verdict"]
mkpath(DST)
cp(joinpath(SRC, "index.html"), joinpath(DST, "index.html"); force = true)
for st in STAGES
    s = joinpath(SRC, st)
    isfile(joinpath(s, "report.html")) || (println("skip $st (no report)"); continue)
    d = mkpath(joinpath(DST, st))
    for f in readdir(d)                                # drop stale figures
        (endswith(f, ".png") || (endswith(f, ".html") && f != "report.html")) && rm(joinpath(d, f))
    end
    for f in readdir(s)
        (endswith(f, ".html") || endswith(f, ".png")) && cp(joinpath(s, f), joinpath(d, f); force = true)
    end
    println("published $st")
end
