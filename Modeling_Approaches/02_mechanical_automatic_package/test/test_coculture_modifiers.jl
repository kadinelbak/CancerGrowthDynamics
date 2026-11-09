# Coculture modifiers of the zoo (zoo_models.jl L.cc, Stage 3a, 2026-11-08): kill factor φ, protection ρ,
# cytostatic γ (immediate or through the drug delay chain) and interference η. Checks the zoo right-hand
# side against the written-out equations, and that a line without modifiers (or with neutral ones) is unchanged.
#   julia --project=. test/test_coculture_modifiers.jl
using Test
isdefined(Main, :U_R) || (const U_R = (10 - 0.5) / (10 - 0.05))
isdefined(Main, :EC50) || (const EC50 = 1.0)
isdefined(Main, :line_spec) || include(joinpath(@__DIR__, "..", "zoo_models.jl"))

# two logistic lines with a 3-compartment linear kill and the Norton–Simon factor, as the carried Stage 2b models
const PR = (rN=0.83, rC=0.80, K=5300.0, βN=0.80, kN=0.92, βC=0.49, kC=0.50, αNC=0.9, αCN=1.15, dN=0.2, dC=0.3)
function cc_sys(ccN, ccC; D=1.0)
    fx(r, β, k, d) = Dict(:r => r, :K => PR.K, :β => β, :k => k, :d => d, :EC50 => 1.0)
    LN = merge(line_spec(:logistic, :chain3, :loss, :D; fixed=fx(PR.rN, PR.βN, PR.kN, PR.dN)), (ns=true,))
    LC = merge(line_spec(:logistic, :chain3, :loss, :D; fixed=fx(PR.rC, PR.βC, PR.kC, PR.dC)), (ns=true,))
    ccN === nothing || (LN = merge(LN, (cc=ccN,)))
    ccC === nothing || (LC = merge(LC, (cc=ccC,)))
    alpha = [[Slot(0, 1.0), Slot(0, PR.αNC)], [Slot(0, PR.αCN), Slot(0, 1.0)]]
    build_system([LN, LC], [(D=D, R=0.0, lines=[1, 2])], alpha)
end
function rhs(sys, u, p=Float64[])
    du = zeros(length(u)); make_system_model(sys)(du, u, p, 0.0); du
end
# state: [N, C, A_N1..A_N3, A_C1..A_C3]
const U = [1200.0, 800.0, 0.3, 0.5, 0.6, 0.2, 0.4, 0.45]

@testset "coculture modifiers" begin
    base = rhs(cc_sys(nothing, nothing), U)
    @test rhs(cc_sys(CC_NEUTRAL, CC_NEUTRAL), U) ≈ base atol = 1e-12

    # all modifiers on, free slots for the sensitive line and fixed ones for the resistant line
    p = [0.6, 2.0, 0.3, 0.05]                     # φ_N, ρ_N, γ_N, η_N
    ccN = (kf=Slot(1, 0.0), pr=Slot(2, 0.0), cs=Slot(3, 0.0), csimm=false, di=Slot(4, 0.0))
    ccC = (kf=Slot(0, 0.4), pr=Slot(0, 0.0), cs=Slot(0, 0.2), csimm=true, di=Slot(0, 0.1))
    for D in (0.0, 1.0)
        du = rhs(cc_sys(ccN, ccC; D=D), U, p)
        N, C = U[1], U[2]; EN = U[5]; EC = U[8]
        qN = (N + PR.αNC * C) / PR.K; qC = (C + PR.αCN * N) / PR.K
        fkN = p[1] / (1 + p[2] * C / PR.K); fkC = 0.4
        dN = max(1 - p[3] * EN, 0) * PR.rN * N * (1 - qN) - (PR.dN + p[4] * C / PR.K + fkN * PR.βN * EN * (1 - qN)) * N
        dC = max(1 - 0.2 * D, 0) * PR.rC * C * (1 - qC) - (PR.dC + 0.1 * N / PR.K + fkC * PR.βC * EC * (1 - qC)) * C
        @test du[1] ≈ dN rtol = 1e-12
        @test du[2] ≈ dC rtol = 1e-12
        @test du[3:end] ≈ rhs(cc_sys(nothing, nothing; D=D), U)[3:end]      # the delay chains are untouched
    end

    # a cytostatic effect through the delay keeps the chain in the equations when the kill is off
    LN0 = merge(line_spec(:logistic, :chain3, :loss, :D; fixed=Dict(:r => 0.8, :K => PR.K, :β => 0.0, :k => 0.9, :EC50 => 1.0)),
                (cc=(kf=Slot(0, 1.0), pr=Slot(0, 0.0), cs=Slot(1, 0.0), csimm=false, di=Slot(0, 0.0)),))
    tex = join(line_tex(LN0, "N", "N+\\alpha_{NC}C", "K"), "\n")
    @test occursin("\\gamma_{N}^{cs}", tex) && occursin("A_{N,3}", tex) && occursin("dA_{N,1}", tex)

    # freezing a fitted line freezes its modifier slots too
    fz = freeze_cc(ccN, p)
    @test fz.kf == Slot(0, 0.6) && fz.pr == Slot(0, 2.0) && fz.cs == Slot(0, 0.3) && fz.di == Slot(0, 0.05) && fz.csimm == false
end

# phenotype switching (L.sw): sensitive cells switch into Z, counted in N, growing at r_C and killed by β_C
@testset "phenotype switching" begin
    fx(r, β, k, d) = Dict(:r => r, :K => PR.K, :β => β, :k => k, :d => d, :EC50 => 1.0)
    LN = merge(line_spec(:logistic, :chain3, :loss, :D; fixed=fx(PR.rN, PR.βN, PR.kN, PR.dN)), (ns=true,))
    LC = merge(line_spec(:logistic, :chain3, :loss, :D; fixed=fx(PR.rC, PR.βC, PR.kC, PR.dC)), (ns=true,))
    alpha = [[Slot(0, 1.0), Slot(0, PR.αNC)], [Slot(0, PR.αCN), Slot(0, 1.0)]]
    sw0 = (μ0=Slot(0, 0.0), μD=Slot(0, 0.0), r=LC.s.r, β=LC.s.β)
    LNs = merge(LN, (sw=merge(sw0, (μ0=Slot(1, 0.0), μD=Slot(2, 0.0))),))
    sys0 = build_system([LN, LC], [(D=1.0, R=0.0, lines=[1, 2])], alpha)
    sys = build_system([LNs, LC], [(D=1.0, R=0.0, lines=[1, 2])], alpha)
    @test sys.n == sys0.n + 1
    # N's block: [A_N1..A_N3, Z]; C's block follows
    u = [1200.0, 800.0, 0.3, 0.5, 0.6, 150.0, 0.2, 0.4, 0.45]
    du = rhs(sys, u, [0.02, 0.1])
    # no switching and no switched cells: identical to the line without sw
    u0 = copy(u); u0[6] = 0.0
    @test rhs(sys, u0, [0.0, 0.0])[[1, 2, 3, 4, 5, 7, 8, 9]] ≈ rhs(sys0, u0[[1, 2, 3, 4, 5, 7, 8, 9]]) rtol = 1e-12
    N, C, Z = u[1], u[2], u[6]; S = N - Z; EN = u[5]
    qN = (N + PR.αNC * C) / PR.K; μ = 0.02 + 0.1 * EN
    dS = PR.rN * S * (1 - qN) - (PR.dN + PR.βN * EN * (1 - qN) + μ) * S
    dZ = PR.rC * Z * (1 - qN) - (PR.dN + PR.βC * EN * (1 - qN)) * Z + μ * S
    @test du[6] ≈ dZ rtol = 1e-12
    @test du[1] ≈ dS + dZ rtol = 1e-12
    tex = join(line_tex(LNs, "N", "N+\\alpha_{NC}C", "K"), "\n")
    @test occursin("Z_{N}", tex) && occursin("\\mu_D", tex) && occursin("\\beta_{C}", tex)
    fz = freeze_sw(LNs.sw, [0.02, 0.1])
    @test fz.μ0 == Slot(0, 0.02) && fz.μD == Slot(0, 0.1)
end
