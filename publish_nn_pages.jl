# Copies the neural-network / UDE track's report pages (HTML + images) into docs/ so the homepage's
# "Older reports" footer can link to them on GitHub Pages. Rerun after the NN track regenerates them:
#   julia publish_nn_pages.jl

const SRC = joinpath(@__DIR__, "Modeling_Approaches", "03_automatic_neural_network", "outputs")
const DST = joinpath(@__DIR__, "docs", "outputs", "nn_track")

n = 0
for (dir, _, files) in walkdir(SRC), f in files
    any(endswith(f, e) for e in (".html", ".png", ".svg")) || continue
    rel = relpath(joinpath(dir, f), SRC)
    mkpath(dirname(joinpath(DST, rel)))
    cp(joinpath(dir, f), joinpath(DST, rel); force=true)
    global n += 1
end
println("Copied $n files to $DST")
