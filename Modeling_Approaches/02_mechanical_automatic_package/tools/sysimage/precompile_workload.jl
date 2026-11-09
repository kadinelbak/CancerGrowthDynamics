# A short run of what the pipeline steps do, so the image holds their compiled code (build_sysimage.jl).
using OrdinaryDiffEq, CSV, DataFrames, Plots, ForwardDiff
f!(du, u, p, t) = (du[1] = p[1] * u[1] * (1 - u[1] / p[2]); nothing)
sol = solve(ODEProblem(f!, [100.0], (0.0, 14.0), [0.8, 4000.0]), Tsit5(); reltol=1e-8, abstol=1e-8)
ForwardDiff.gradient(p -> solve(ODEProblem(f!, [100.0], (0.0, 14.0), p), Tsit5())(14.0)[1], [0.8, 4000.0])
df = DataFrame(t=sol.t, y=first.(sol.u))
io = IOBuffer(); CSV.write(io, df); CSV.read(IOBuffer(take!(io)), DataFrame)
plt = plot(df.t, df.y; label="logistic"); scatter!(plt, df.t, df.y); savefig(plt, tempname() * ".png")
