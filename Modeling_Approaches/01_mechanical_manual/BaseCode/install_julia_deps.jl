using Pkg

println("Installing Julia dependencies for notebooks and optimization...")

# Core notebook support
Pkg.add("IJulia")

# Optimization stack
# - Optimization.jl: common interface
# - OptimizationOptimJL.jl: bridge to Optim.jl-based algorithms
# - Optim.jl: NelderMead, BFGS, etc.
Pkg.add(["Optimization", "OptimizationOptimJL", "Optim"])

# Optional but often useful for AD with Optimization.jl
try
	Pkg.add("ForwardDiff")
catch e
	@warn "ForwardDiff install failed (optional)" exception=e
end

println("Precompiling packages (this may take a moment)...")
Pkg.precompile()

println("All Julia dependencies installed and precompiled.")