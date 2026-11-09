# Precompiled Julia image for the pipeline (user, 2026-11-09): every pipeline step starts a fresh Julia process,
# and each spent ~30–40 s loading and compiling the packages (a report step took ~41 s in all). With this image a
# step starts in a few seconds. It changes no result: the same packages, only compiled ahead of time.
#
#   julia --project=tools/sysimage tools/sysimage/build_sysimage.jl        (once, ~10–20 min; again after a
#                                                                            package update, run_pipeline.jl says so)
#
# Writes tools/sysimage/a2780_sysimage.<so|dll|dylib> and manifest.sha256 (the package Manifest it was built
# from). run_pipeline.jl uses the image only while that Manifest is unchanged; otherwise it runs without it.
using Pkg
Pkg.instantiate()
using PackageCompiler, SHA

const PKG = normpath(joinpath(@__DIR__, "..", ".."))
const EXT = Sys.iswindows() ? "dll" : Sys.isapple() ? "dylib" : "so"
const IMG = joinpath(@__DIR__, "a2780_sysimage.$EXT")
# the heavy packages every step loads (Project.toml [deps]); the small standard libraries load fast anyway
const PKGS = [:GrowthParameterEstimation, :OrdinaryDiffEq, :ForwardDiff, :BlackBoxOptim, :CSV, :DataFrames,
              :XLSX, :JSON3, :Plots, :HTTP]
create_sysimage(PKGS; sysimage_path=IMG, project=PKG,
                precompile_execution_file=joinpath(@__DIR__, "precompile_workload.jl"))
write(joinpath(@__DIR__, "manifest.sha256"), bytes2hex(open(sha256, joinpath(PKG, "Manifest.toml"))))
println("Wrote $IMG")
