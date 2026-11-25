# Control plant for optimal control of cisplatin dosing, built from the selected compartment models.
#
# The biology is NOT refitted here. The cell equations are the zoo's own right-hand side
# (zoo_models.jl: line_rhs!; simulated through a flattened, type-stable copy that the tests check against
# it to 1e-12) with the parameters of the carried-forward models in outputs/zoo_fits/zoo_state.jls.
# Only the dose input and a toxicity proxy are added:
#
#   state x = [N, C, A_N,1..3, A_C,1..3, (low resource: S_C,1..3), B, T]
#     N, C      sensitive (A2780) and resistant (A2780cis) cells (Mean Cells units)
#     A_X,j     drug-effect compartments of line X (3-compartment chain, rate k_X, driven by u)
#     S_C,j     low-resource stress compartments of A2780cis (low-resource scenario only)
#     B         systemic drug burden, dB/dt = s_B·u − d_B·B               (de Pillis & Radunskaya 2001)
#     T         host toxicity ([toxicity] model): xie2020 (default) = kidney excess E, dE/dt = k_N·Emax·γu/(EC50 + γu) − k_N·E
#               (Xie et al. 2020); depillis = toxicity burden 1 − H, dT/dt = a_H(1 − e^(−B))(1 − T) − r_H·T(1 − T)
#   control u = cisplatin in the medium (µM), piecewise constant over each decision interval
#
# Scenarios for the drug response in coculture (user decision, 2026-09-29):
#   :A  the Stage 3 optimal model as selected (coculture fit; the resistant kill is inherited from
#       Stage 2, refitted, or fixed at 0 when the selected candidate has no resistant cisplatin term:
#       stage3_betaC() says which)
#   :B  Stage 2 monoculture drug response for both lines, crowding fitted to untreated cocultures
#       only (the variant that predicted the held-out treated cocultures best)

include(joinpath(@__DIR__, "zoo_common.jl"))
include(joinpath(@__DIR__, "zoo_models.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization, TOML, ForwardDiff, Statistics
import LinearAlgebra
const LA = LinearAlgebra

const CFG = TOML.parsefile(joinpath(@__DIR__, "control_config.toml"))
# ZOO_STATE / AN_STATE override the fitted states (e.g. a backed-up earlier fit, to compare control runs)
const ZOO_STATE = deserialize(get(ENV, "ZOO_STATE", joinpath(@__DIR__, "outputs", "zoo_fits", "zoo_state.jls")))
const AN_STATE = deserialize(get(ENV, "AN_STATE", joinpath(@__DIR__, "outputs", "zoo_analyses", "analyses_state.jls")))
zgroup(id) = ZOO_STATE.groups[findfirst(g -> g.id == id, ZOO_STATE.groups)]

# ── reuse of finished control runs (2026-11-01) ──────────────────────────────────────────────────────────
# A control step (run_control.jl, run_control_fasting.jl, run_control_sweep.jl) is skipped when everything it
# reads is unchanged since its last finished run: the control code and configs, the data, the package
# versions, and the parts of the fitted states it uses (each group's carried model and its bootstrap draws,
# the Ce0 crowding fit). A refit that leaves the carried models and draws the same therefore costs no control
# time; anything else reruns. CONTROL_FORCE=1 (run_pipeline.jl --recompute-control) always recomputes.
# The fingerprint is a SHA-256 written next to the outputs after a successful run.
const CONTROL_CODE = ["control_model.jl", "control_mpc.jl", "control_fasting.jl", "zoo_common.jl", "zoo_models.jl",
                      "zoo_data.jl", "bootstrap_tools.jl", "stage_names.jl", "control_config.toml", "stages.toml",
                      "Manifest.toml", joinpath("data", "schedules.toml"), joinpath("data", "measurements.csv")]
function control_fingerprint(script::AbstractString, extra::Pair...)
    io = IOBuffer()
    for f in vcat(CONTROL_CODE, basename(script))
        fp = joinpath(@__DIR__, f)
        write(io, f, 0x00, isfile(fp) ? read(fp) : UInt8[], 0x00)
    end
    for g in ZOO_STATE.groups
        c = g.cands[1]
        write(io, string(g.id), 0x00, string(c.name), 0x00, join(string.(c.pnames), ","), 0x00, Float64.(c.params))
        if hasproperty(c, :fixedcc)     # sorted: a Dict's order is not stable between sessions
            fc = c.fixedcc isa AbstractDict ? sort!([string(k) => v for (k, v) in c.fixedcc]; by=first) : c.fixedcc
            write(io, repr(fc), 0x00)
        end
        if haskey(ZOO_STATE.results, g.id)
            P = ZOO_STATE.results[g.id].res.pars[1]
            write(io, Int64(size(P, 1)), Int64(size(P, 2)), Float64.(vec(P)))
        end
    end
    write(io, repr(ZOO_STATE.B), 0x00, Float64.(AN_STATE[:validation].ce1[1].params))
    for (k, v) in extra
        write(io, string(k), "=", string(v), 0x00)
    end
    bytes2hex(sha256(take!(io)))
end
# true when `fpfile` holds this fingerprint and every output in `outs` exists (and CONTROL_FORCE is not 1)
control_unchanged(fpfile, fp, outs...) = get(ENV, "CONTROL_FORCE", "0") != "1" && isfile(fpfile) &&
    strip(read(fpfile, String)) == fp && all(isfile, outs)

# resistant kill β_C of the carried Stage 3 model, and where it comes from: refitted, dropped (the
# "no resistant cisplatin term" candidates fix it at 0) or inherited from Stage 2b
function stage3_betaC()
    c3 = zgroup("s3").cands[1]; s = c3.sys.lines[2].s.β
    dropped = s.i == 0 && hasproperty(c3, :fixedcc) && haskey(c3.fixedcc, :βC)
    s.i > 0 ? (value=c3.params[s.i], free=true, dropped=false,
               source="Stage 3 (refitted" * (c3.bounds[s.i][1] > 0 ? ", floor $(round(c3.bounds[s.i][1]; sigdigits=3)))" : ", no floor)")) :
    dropped ? (value=s.v, free=false, dropped=true, source="Stage 3 (no resistant cisplatin term in coculture: selected, β_C = 0)") :
              (value=s.v, free=false, dropped=false, source="Stage 2 (inherited)")
end

# ── flattened cell equations ─────────────────────────────────────────────────
# Flattened, type-stable copy of the plant's cell equations (the same modules and values as the zoo's
# line_rhs! / line_rhs_multi!, which are generic and dispatch dynamically, about 50× slower). The tests
# check it against the zoo right-hand side to 1e-12 (normal medium: test_control.jl; with the low-resource
# modules: test_control_fasting.jl). Only the module types of the carried models are supported;
# line_specs stops with an error on anything else rather than silently differing.
struct ModSpec
    input_D::Bool; kernel_chain::Int; k::Float64; o::Int     # chain length (comp1 = 1, chain2 = 2, chain3 = 3; 0 = none/auc/second), rate, first hidden index − 1
    action::Symbol; β::Float64; h::Float64; c::Float64; d0::Float64; ec50::Float64
    ε::Float64; εj::Int                                       # fixed strength, or index into the free strengths (0 = fixed)
    kernel::Symbol; ω::Float64; ζ::Float64                    # delay kernel (zoo_models.jl KERNELS) and second-order ω, ζ
end
struct LineSpec
    ix::Int; other::Int; α::Float64; logistic::Bool; r::Float64; K::Float64; θ::Float64; d::Float64
    mods::Vector{ModSpec}
    ns::Bool                    # Norton–Simon: drug losses × growth factor (zoo_models.jl growth_factor)
    chi::Float64                # drug × low-resource interaction χ: drug losses × (1 + χ·u_R/U_R) (zoo line_rhs_multi!)
    # coculture modifiers (zoo_models.jl L.cc; neutral 1, 0, 0, true, 0 for lines without them): kill factor φ,
    # protection ρ, cytostatic γ (immediate or through the first drug module's delay), interference η
    kf::Float64; pr::Float64; cs::Float64; csimm::Bool; di::Float64
    pn::Float64                 # Norton–Simon exponent (zoo_models.jl nspow; 1 = full coupling)
end
const FAST_KERNELS = (:none, :comp1, :chain2, :chain3, :second, :auc)
# iε: positions in p of free low-resource strengths (control_fasting.jl); none in normal medium
function line_specs(S, p, iε)
    out = LineSpec[]
    for (j, l) in enumerate(S.envs[1].lines)
        L = S.lines[l]; xi = S.xi[1][j]; off = S.hb[1][j]
        L.growth in (:logistic, :theta) || error("growth law $(L.growth) not supported in the fast plant")
        sw_on(L) && error("phenotype switching (L.sw) not supported in the fast plant")
        mods = ModSpec[]
        for M in modules(L)
            ncell(M.action) == 0 && naux(M.action) == 0 || error("action $(M.action) not supported in the fast plant")
            M.kernel in FAST_KERNELS || error("kernel $(M.kernel) not supported in the fast plant")
            M.action in (:hill, :capacity, :loss) || error("action $(M.action) not supported in the fast plant")
            nk = nker(M.kernel); gs(k) = hasproperty(M.s, k) ? sv(getfield(M.s, k), p, Float64) : 0.0
            εs = M.s.ε; εj = εs.i > 0 ? findfirst(==(εs.i), iε) : 0
            εs.i > 0 && εj === nothing && error("unexpected free slot")
            ch = M.kernel in (:comp1, :chain2, :chain3) ? nk : 0
            push!(mods, ModSpec(M.input == :D, ch, ch > 0 ? gs(:k) : 0.0, off, M.action, gs(:β), gs(:h), gs(:c), gs(:d0), ec50_of(M.s, p, Float64),
                                εs.i > 0 ? NaN : εs.v, εj, M.kernel, gs(:ω), gs(:ζ)))
            off += mwidth(M)
        end
        m = 3 - j
        cc = cc_on(L) ? L.cc : CC_NEUTRAL
        (cc_on(L) && L.input != :D) && error("coculture modifiers on a line whose own module is not the drug module")
        push!(out, LineSpec(xi + 1, S.xi[1][m] + 1, sv(S.alpha[l][S.envs[1].lines[m]], p, Float64), L.growth == :logistic,
                            sv(L.s.r, p, Float64), sv(L.s.K, p, Float64), sv(L.s.θ, p, Float64), sv(L.s.d, p, Float64), mods, ns_on(L),
                            hasproperty(L, :chi) ? sv(L.chi, p, Float64) : 0.0,
                            sv(cc.kf, p, Float64), sv(cc.pr, p, Float64), sv(cc.cs, p, Float64), cc.csimm, sv(cc.di, p, Float64),
                            pn_of(L.s, p, Float64)))
    end
    out
end
# Cell (and drug/stress chain) derivatives for dose u and low-resource level uR. f = ε factor of the dose,
# εv = free low-resource strengths, koff/v = refeeding fade rate (NaN = fitted rate) and fasting input.
function spec_cells_rhs!(du, x, u, uR, f, spec, εv, koff, v)
    T = promote_type(eltype(x), typeof(u), typeof(uR)); z = zero(T)
    @inbounds for L in spec
        X = max(x[L.ix], z); Y = max(x[L.other], z); crowd = X + L.α * Y
        Keff = T(L.K); loss = T(L.d); dloss = z                   # dloss: drug-module losses (Norton–Simon scales these)
        Ecs = z; first_drug = true                                # cytostatic exposure: the first drug module's dose or effect state
        for M in L.mods
            ε = M.εj == 0 ? T(M.ε) : 1 + (εv[M.εj] - 1) * f
            wi = (M.input_D ? T(u) : T(uR)) * ε
            E = wi; o = M.o
            if M.kernel_chain > 0                    # effect compartment / transit chain of n compartments at rate k
                # low-resource stress may fade faster after refeeding (sensitivity only; default = fitted rate)
                k = (M.input_D || isnan(koff)) ? T(M.k) : M.k * v + koff * (1 - v)
                du[o+1] = k * (wi - x[o+1])
                for i in 2:M.kernel_chain
                    du[o+i] = k * (x[o+i-1] - x[o+i])
                end
                E = max(x[o+M.kernel_chain], z)
            elseif M.kernel == :auc                  # cumulative exposure: never fades in the fitted model; the
                fade = (M.input_D || isnan(koff)) ? z : koff * (1 - v)       # refeeding sensitivity adds a fade when fed
                du[o+1] = wi - fade * x[o+1]
                E = max(x[o+1], z)
            elseif M.kernel == :second               # second-order (damped oscillator) delay
                du[o+1] = x[o+2]; du[o+2] = M.ω^2 * (wi - x[o+1]) - 2M.ζ * M.ω * x[o+2]
                E = max(x[o+1], z)
            end
            if M.input_D && first_drug
                Ecs = L.csimm ? wi : E; first_drug = false
            end
            if M.action == :hill
                lo = E > 0 ? M.β * E^M.h / (T(M.ec50)^M.h + E^M.h) : z
                M.input_D ? (dloss += lo) : (loss += lo)
            elseif M.action == :capacity
                Keff /= 1 + M.c * E
            else                                     # :loss
                lo = M.d0 * wi + M.β * E
                M.input_D ? (dloss += lo) : (loss += lo)
            end
        end
        q = max(crowd / Keff, z)
        y = Y / T(L.K)                              # coculture modifiers (neutral unless Stage 3a carries them)
        dloss *= L.kf / (1 + L.pr * y)
        loss += L.di * y
        gm = L.cs == 0 ? one(T) : max(1 - L.cs * Ecs, z)
        dloss *= 1 + L.chi * T(uR) / T(U_R)        # drug × low-resource interaction (0 unless Stage 4 carries χ)
        loss += L.ns ? dloss * nspow(L.logistic ? max(1 - q, z) : max(1 - q^L.θ, z), L.pn) : dloss
        du[L.ix] = gm * L.r * X * (L.logistic ? 1 - q : 1 - q^L.θ) - loss * X
    end
    nothing
end
const NO_EPS = Float64[]

# ── plant ─────────────────────────────────────────────────────────────────────

struct Plant
    name::String
    sys::Any                  # zoo system (one environment, lines [N, C])
    p::Vector{Float64}        # free-slot parameter vector of the zoo lines
    nx::Int
    iN::Int; iC::Int; iB::Int; iT::Int
    K::Float64                # shared carrying capacity (for normalisation and the cancer bound)
    tox::NamedTuple           # default_tox(): r_H, a_H, d_B, s_B (depillis) and xie, k_N, Emax, EC50, γ, Elim (xie2020)
    uR::Float64               # low-resource level of the environment (0 in normal medium)
    fp::Any                   # low-resource environment: the two-input plant of control_fasting.jl held at v = 1, else nothing
    spec::Any                 # normal medium: flattened cell equations (spec_cells_rhs!), or nothing = the zoo right-hand side
end
line_specs(base::Plant, iε) = line_specs(base.sys, base.p, iε)

# replace fixed slot values of a line (e.g. a bootstrap draw of carried parameters)
function set_slots(L, d::Dict{Symbol,Float64})
    s = L.s
    for (k, v) in d
        s = merge(s, NamedTuple{(k,)}((Slot(0, v),)))
    end
    merge(L, (s=s,))
end

# ── host toxicity model (control_config.toml [toxicity] model) ─────────────────────────────────────────────
# "depillis": healthy-tissue loss T driven by the drug burden B (de Pillis & Radunskaya 2001), limit B_tox.
# "xie2020" (default since 2026-11-21, user): kidney toxicity, the indirect-response nephrotoxicity model of Xie et
# al. 2020 (AAPS J 22:96) fitted in 50 ovarian cancer patients on cisplatin, dResp/dt = kin(1 + Emax·Cp/(EC50 + Cp))
# − kout·Resp with kin = kout. The state T holds the excess E = Resp − 1 (0 at baseline):
#   dE/dt = k_N·Emax·Cp/(EC50 + Cp) − k_N·E,   Cp = γ·u
# Cp is a surrogate plasma level: the medium concentration used as a stand-in for plasma level, scaled by γ. γ is
# set at load time so that constant 1 µM gives the same time-averaged excess Ē as the clinical schedule (75 mg/m²
# every 21 days, a 1-hour infusion; plasma curve of intact cisplatin scaled from Liston & Davis 2017, Clin Cancer
# Res 23:3489); the limit E_lim is the peak E of that schedule at its periodic steady state. B stays in the state
# vector but does not feed the toxicity under xie2020.
const TOX_MODEL = String(get(CFG["toxicity"], "model", "depillis"))
TOX_MODEL in ("depillis", "xie2020") || error("control_config.toml [toxicity] model must be \"depillis\" or \"xie2020\", got $TOX_MODEL")

# the clinical calibration of the xie2020 model (hours inside, k_N in /day): ΔE per cycle, mean Ē and peak E at the
# periodic steady state of the clinical schedule, and the γ that gives Ē under a constant 1 µM
function xie_calibration(t=CFG["toxicity"])
    kN = Float64(t["k_N"]); Emax = Float64(t["Emax"]); EC50 = Float64(t["EC50_uM"])
    Cmax = Float64(t["clin_Cmax_uM"]) * Float64(t["clin_dose_mg_m2"]) / Float64(t["clin_ref_dose_mg_m2"])
    ke = log(2) / Float64(t["clin_half_life_h"]); tinf = Float64(t["clin_infusion_h"]); cyc = Float64(t["clin_cycle_days"])
    R = Cmax * ke / (1 - exp(-ke * tinf))                     # zero-order infusion rate giving Cmax at the end of the infusion
    Cp(th) = th <= tinf ? R / ke * (1 - exp(-ke * th)) : Cmax * exp(-ke * (th - tinf))
    k = kN / 24; dt = 0.005; nstep = round(Int, cyc * 24 / dt)
    f(E, th) = k * (Emax * Cp(th) / (EC50 + Cp(th)) - E)
    E = 0.0; ΔE = NaN; Ebar = NaN; Epk = NaN
    for c in 1:40                                             # 40 cycles: e^(−k_N·21·40) ≈ 0, a periodic steady state
        s = 0.0; pk = E
        for i in 0:nstep-1                                    # RK4 on the hour scale
            th = i * dt
            k1 = f(E, th); k2 = f(E + dt / 2 * k1, th + dt / 2); k3 = f(E + dt / 2 * k2, th + dt / 2); k4 = f(E + dt * k3, th + dt)
            E += dt / 6 * (k1 + 2k2 + 2k3 + k4); s += E * dt; pk = max(pk, E)
        end
        c == 1 && (ΔE = pk)                                   # first cycle from E = 0: the rise per dose
        Ebar = s / (nstep * dt); Epk = pk
    end
    γ = Ebar * EC50 / (Emax - Ebar)                           # Emax·γ/(EC50 + γ) = Ē
    (kN=kN, Emax=Emax, EC50=EC50, γ=γ, ΔE=ΔE, Ebar=Ebar, Elim=Epk)
end
const XIE = TOX_MODEL == "xie2020" ? xie_calibration() : nothing
# CONTROL_TAG (2026-11-21): a side run writes to outputs/control<TAG>/ (14 days) or outputs/control<TAG>_<days>d/, and its
# reports skip the docs/ copies and the site summaries, so the main control results stay as they are
const CONTROL_TAG = get(ENV, "CONTROL_TAG", "")
control_subdir(days, base_days) = days == base_days ? "control$(CONTROL_TAG)" : "control$(CONTROL_TAG)_$(round(Int, days))d"
XIE === nothing || println("Kidney toxicity (Xie et al. 2020) calibrated to 75 mg/m² every 21 days: γ = $(round(XIE.γ; sigdigits=4)) per µM, ΔE = $(round(XIE.ΔE; sigdigits=4)), Ē = $(round(XIE.Ebar; sigdigits=4)), E_lim = $(round(XIE.Elim; sigdigits=4))")

# γ_mult: the γ multiplier of the robustness sweep (xie2020 only)
function default_tox(; γ_mult=1.0)
    t = CFG["toxicity"]
    base = (r_H=Float64(t["r_H"]), a_H=Float64(t["a_H"]), d_B=Float64(t["d_B"]), s_B=Float64(t["s_B"]))
    XIE === nothing && return merge(base, (xie=false, kN=0.0, Emax=0.0, EC50=1.0, γ=0.0, Elim=Float64(t["B_tox"])))
    merge(base, (xie=true, kN=XIE.kN, Emax=XIE.Emax, EC50=XIE.EC50, γ=XIE.γ * γ_mult, Elim=XIE.Elim))
end
# the toxicity limit the objective, barriers and violation flags use: E_lim (xie2020) or B_tox (depillis)
tox_limit(tox=default_tox()) = tox.Elim
# the kidney excess drive k_N·Emax·Cp/(EC50 + Cp) at dose u (µM in the medium)
xie_drive(tp, u) = (Cp = tp.γ * max(u, zero(u)); tp.kN * tp.Emax * Cp / (tp.EC50 + Cp))

# report texts for the active toxicity model (control_report.jl, control_fasting_report.jl)
# LaTeX of the objective with the active toxicity state and limit (T/B_tox → E/E_lim under xie2020)
tox_tex(s) = TOX_MODEL == "xie2020" ? replace(s, "{T}{B_{tox}}" => "{E}{E_{lim}}", "{T}{B_{tox}" => "{E}{E_{lim}") : s
const XIE_REF = "Xie F, et al. (2020). Population pharmacokinetics and nephrotoxicity of cisplatin in ovarian cancer patients (indirect-response model of serum creatinine). AAPS J 22:96, doi 10.1208/s12248-020-00489-2."
const LISTON_REF = "Liston DR, Davis M (2017). Clinically relevant concentrations of anticancer drugs: a guide for nonclinical studies. Clin Cancer Res 23:3489–3498."
const XIE_CALIB = XIE === nothing ? "" :
    "Calibration (ASSUMPTION): γ = $(round(XIE.γ; sigdigits=3)) per µM, so that a constant 1 µM gives the time-averaged kidney excess Ē = $(round(XIE.Ebar; sigdigits=3)) of the clinical " *
    "schedule (75 mg/m² every 21 days as a 1-hour infusion, intact-cisplatin plasma curve scaled from Liston & Davis 2017); one clinical dose raises E by ΔE = $(round(XIE.ΔE; sigdigits=3)), " *
    "and the limit E_lim = $(round(XIE.Elim; sigdigits=3)) is the peak E of that schedule at its periodic steady state."
const TX = TOX_MODEL == "xie2020" ? (
    name="Kidney toxicity E", short="kidney excess E", sym="E", lim="E_lim", limtex="E_{lim}", unit="excess over baseline", ylabel="E (excess over baseline)",
    eq=raw"\frac{dE}{dt}=k_N\,E_{max}\frac{\gamma u}{EC_{50}+\gamma u}-k_N\,E",
    eqf=raw"\frac{dE}{dt}=k_N\,E_{max}\frac{\gamma u}{EC_{50}+\gamma u}\left(1-p_H\,P\right)-k_N\,E",
    what="E = the excess of the kidney-injury response over its baseline (0 = none), from the cisplatin nephrotoxicity model of Xie et al. 2020, fitted in 50 ovarian cancer patients. " *
         "The medium concentration u is used as a stand-in for plasma level, Cp = γ·u. " * XIE_CALIB,
    limit_note="The dashed line is the limit E_lim = $(XIE === nothing ? "" : round(XIE.Elim; sigdigits=3)), the peak kidney excess of the standard clinical schedule.",
    caveat="Nothing in these in-vitro data measures host toxicity. The kidney model and its values come from patients (Xie et al. 2020); the medium concentration is used as a stand-in for plasma level, " *
           "and the scale γ and the limit E_lim are calibrated to the clinical schedule (an assumption, swept ×1, ×3, ×10). Treat toxicity-limited results as illustrative.",
    refs=[XIE_REF, LISTON_REF]) : (
    name="Toxicity burden T", short="toxicity burden T", sym="T", lim="B_tox", limtex="B_{tox}", unit="fraction 0–1", ylabel="T (fraction)",
    eq=raw"\frac{dB}{dt}=s_B\,u-d_B\,B,\qquad \frac{dT}{dt}=a_H\left(1-e^{-B}\right)(1-T)-r_H\,T(1-T)",
    eqf=raw"\frac{dT}{dt}=a_H\,\left(1-p_H\,P\right)\left(1-e^{-B}\right)(1-T)-r_H\,T\,(1-T)",
    what="T = 1 − H = fraction of healthy host tissue lost (literature model of de Pillis & Radunskaya 2001), driven by the systemic drug burden B.",
    limit_note="The dashed line is the bound B_tox ($(CFG["toxicity"]["B_tox"]): healthy tissue kept above $(round(Int, 100(1 - CFG["toxicity"]["B_tox"])))%).",
    caveat="Nothing in these in-vitro data measures host toxicity. Its structure and values come from de Pillis & Radunskaya (2001); the mapping from medium concentration to systemic burden (s_B) is an assumption. Treat toxicity-limited results as illustrative.",
    refs=String[])

# fast = true: simulate the cells with the flattened equations (normal medium only; a model with modules the
# flattened copy does not support falls back to the zoo right-hand side, which gives the same values, slower)
function make_plant(name, lines, alpha, p; uR=0.0, tox=default_tox(), fast=true)
    sys = build_system(lines, [(D=0.0, R=uR, lines=[1, 2])], alpha)
    n = sys.n
    Kslot = lines[1].s.K
    K = Kslot.i > 0 ? p[Kslot.i] : Kslot.v
    spec = (fast && uR == 0) ? (try line_specs(sys, p, Int[]) catch; nothing end) : nothing
    Plant(name, sys, p, n + 2, sys.xi[1][1] + 1, sys.xi[1][2] + 1, n + 1, n + 2, K, tox, uR, nothing, spec)
end

# The Stage 3 coculture lines and parameters for a scenario (point estimates or a bootstrap draw b)
function scenario_plant(scen::Symbol; b::Union{Nothing,Int}=nothing, tox=default_tox())
    g3 = zgroup("s3"); c3 = g3.cands[1]
    LN, LC = c3.sys.lines
    # a carried parameter of an earlier stage by name ("θ" → θ_N or θ_C), from its fit or bootstrap replicate b;
    # nothing when the carried model has no such parameter (e.g. h for a linear kill)
    function draw(gid, sym)
        c = zgroup(gid).cands[1]; i = findfirst(n -> n in ("$(sym)_N", "$(sym)_C"), c.pnames)
        i === nothing && return nothing
        b === nothing ? c.params[i] : ZOO_STATE.results[gid].res.pars[1][b, i]
    end
    setdraw!(d, sym, gid) = (v = draw(gid, sym); v === nothing || (d[sym] = v); d)
    fixed_N = Dict{Symbol,Float64}(); fixed_C = Dict{Symbol,Float64}()
    if b !== nothing       # carried parameters from the same bootstrap replicate of their own stage
        setdraw!(fixed_N, :r, "s1-a2780"); setdraw!(fixed_N, :θ, "s1-a2780"); setdraw!(fixed_C, :r, "s1-a2780cis"); setdraw!(fixed_C, :θ, "s1-a2780cis")
        for s in (:h, :k)
            setdraw!(fixed_N, s, "s2-a2780"); setdraw!(fixed_C, s, "s2-a2780cis")
        end
    end
    if scen == :A
        p = b === nothing ? copy(c3.params) : ZOO_STATE.results["s3"].res.pars[1][b, :]
        # an inherited resistant kill follows its own Stage 2 bootstrap draw, as h and k do
        b !== nothing && LC.s.β.i == 0 && !stage3_betaC().dropped && setdraw!(fixed_C, :β, "s2-a2780cis")
        name = "A: Stage 3 coculture model as selected"
    else
        # crowding (K, α_NC, α_CN and d_N where the model has it) fitted to the untreated cocultures; these
        # lead the Stage 3 parameter vector (run_heldout.jl checks the same), the kills come from Stage 2.
        # The untreated-coculture crowding fit has no bootstrap of its own: a draw b BORROWS the Stage 3a
        # bootstrap ratio of replicate b to its estimate (pars[b] / estimate), applied to the Ce0 values
        # (an approximation, disclosed in the report)
        ce0 = AN_STATE[:validation].ce1[1].params; n = length(ce0)
        p = b === nothing ? copy(ce0) : ZOO_STATE.results["s3"].res.pars[1][b, 1:n] .* (ce0 ./ c3.params[1:n])
        setdraw!(fixed_N, :β, "s2-a2780"); setdraw!(fixed_C, :β, "s2-a2780cis")
        # monoculture drug response: the Stage 2b delay rate (where Stage 3a refitted it) and no coculture drug
        # modifier (kill factor, protection, cytostatic effect; zoo_heldout.jl mono_drug_line does the same)
        LN.s.k.i > 0 && setdraw!(fixed_N, :k, "s2-a2780"); LC.s.k.i > 0 && setdraw!(fixed_C, :k, "s2-a2780cis")
        nodrugcc(L) = cc_on(L) ? merge(L, (cc=merge(L.cc, (kf=Slot(0, 1.0), pr=Slot(0, 0.0), cs=Slot(0, 0.0))),)) : L
        LN, LC = nodrugcc(LN), nodrugcc(LC)
        name = "B: monoculture drug response for both lines, crowding from untreated cocultures"
    end
    any(isnan, p) && return nothing
    make_plant(name, [set_slots(LN, fixed_N), set_slots(LC, fixed_C)], c3.sys.alpha, p; tox=tox)
end

# Low-resource environment: lr_plant(scen) is defined in control_fasting.jl (the two-input plant held at
# v = 1), so the 14-day and the long (dose + fasting) controllers use the same dose-dependent drug ×
# low-resource model.

# continuous-time right-hand side for a given dose u (µM); ref = true forces the zoo's own right-hand side
# (the reference the flattened equations are tested against)
function plant_rhs!(du, x, u, P::Plant; ref=false)
    T = promote_type(eltype(x), typeof(u))
    P.fp === nothing || return fast_rhs!(du, x, (T(u), one(T)), P.fp)      # low resource: fasting plant at v = 1
    if P.spec !== nothing && !ref                                          # normal medium, flattened equations
        spec_cells_rhs!(du, x, T(u), zero(T), zero(T), P.spec, NO_EPS, NaN, 0)
    else
        zoo_cells_rhs!(du, x, u, P, T)
    end
    B = x[P.iB]; Tx = x[P.iT]; tp = P.tox
    du[P.iB] = tp.s_B * u - tp.d_B * B
    du[P.iT] = tp.xie ? xie_drive(tp, u) - tp.kN * Tx :
                        tp.a_H * (1 - exp(-max(B, zero(T)))) * (1 - Tx) - tp.r_H * Tx * (1 - Tx)
    nothing
end
function zoo_cells_rhs!(du, x, u, P::Plant, T)
    S = P.sys; env = S.envs[1]
    for (j, l) in enumerate(env.lines)
        xi = S.xi[1][j]; L = S.lines[l]
        crowd = max(x[xi+1], zero(T)); other = zero(T)
        for (m, lm) in enumerate(env.lines)
            m == j && continue
            y = max(x[S.xi[1][m]+1], zero(T))
            crowd += sv(S.alpha[l][lm], P.p, T) * y; other += y
        end
        if hasproperty(L, :extra)          # Stage 4 lines: drug module + low-resource module
            line_rhs_multi!(du, x, xi, S.hb[1][j], L, P.p, (D=u, R=env.R), crowd, T; other=other)
        else
            w = L.input == :D ? T(u) : T(env.R)
            line_rhs!(du, x, xi, S.hb[1][j], L, P.p, w, crowd, T; other=other)
        end
    end
    nothing
end
function plant_f(x, u, P::Plant)
    T = promote_type(eltype(x), typeof(u))
    du = zeros(T, length(x)); plant_rhs!(du, x, u, P); du
end

# ── discrete transition F(x, u) over one decision interval (fixed-step RK4) ─────

const RK4_H = Float64(CFG["schedule"]["rk4_step_days"])
function F(x, u, P::Plant, Δ; h=RK4_H)
    n = max(1, ceil(Int, Δ / h - 1e-9)); h = Δ / n
    T = promote_type(eltype(x), typeof(u))
    y = T.(x)
    for _ in 1:n
        k1 = plant_f(y, u, P); k2 = plant_f(y .+ (h / 2) .* k1, u, P)
        k3 = plant_f(y .+ (h / 2) .* k2, u, P); k4 = plant_f(y .+ h .* k3, u, P)
        y = y .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    y
end
# dense trajectory of a schedule (u per interval) for plotting
function trajectory(x0, U, P::Plant, Δ; dt=0.1)
    ts = [0.0]; xs = [copy(x0)]; x = copy(x0); t = 0.0
    m = max(1, round(Int, Δ / dt))
    for u in U, _ in 1:m
        x = F(x, u, P, Δ / m); t += Δ / m
        push!(ts, t); push!(xs, x)
    end
    ts, xs
end

# Local linearisation around (x̄, ū) by forward-mode automatic differentiation:
#   A = ∂F/∂x, B = ∂F/∂u, one ForwardDiff Jacobian of F over z = [x; u] (nx + 1 dual directions).
# With ~9–12 states and one input, forward mode needs nx + 1 passes and reverse mode nx passes,
# so the cost is similar; forward mode differentiates through the RK4 loop directly and needs no tape.
# At a dose bound (u ≤ lo or u ≥ hi) B is the one-sided derivative into the feasible side (second-order
# one-sided difference): at u = 0 with no drug in the system the drug-effect states sit exactly on their
# max(E, 0) clamp, where AD takes the clamped branch and misses the dose sensitivity (2026-11-08).
function linearize(x, u, P::Plant, Δ; lo=-Inf, hi=Inf)
    nx = length(x)
    J = ForwardDiff.jacobian(z -> F(z[1:nx], z[nx+1], P, Δ), vcat(x, u))
    Fx = F(x, u, P, Δ); B = J[:, nx+1]
    s = u <= lo ? 1.0 : u >= hi ? -1.0 : 0.0
    s == 0 || (h = 1e-4 * s; B = (4 .* F(x, u + h, P, Δ) .- F(x, u + 2h, P, Δ) .- 3 .* Fx) ./ (2h))
    Fx, J[:, 1:nx], B
end

# Starting counts of every control run: the NOMINAL seeding of the coculture experiments
# (control_config.toml [initial]: 100 cells in total, split by the sensitive:resistant ratio, so
# "25:75" → N = 25, C = 75). Not a measured day-0 value: the zoo fits treat day 0 as the known seeding
# and the candidates' u0 as their starting state, so the controller starts from the design value too.
function nominal_seeding(ratio::AbstractString)
    s, r = parse.(Float64, split(ratio, ":"))
    tot = Float64(CFG["initial"]["seeding_total"])
    (N0=tot * s / (s + r), C0=tot * r / (s + r))
end

# The Stage 3a coculture series with this ratio, treated (Ce1, D > 0) or untreated (Ce0): its environment and the
# u indices of its N and C (the carried model's lines: 1 = A2780, 2 = A2780cis). Found from the nominal seeding.
function coculture_env(ratio::AbstractString, treated::Bool)
    c3 = zgroup("s3").cands[1]; S = c3.sys; obs = sys_obs(S); nom = nominal_seeding(ratio)
    for (e, env) in enumerate(S.envs)
        (input_at(env.D, 0.0) > 0) == treated || continue
        iN, iC = S.xi[e][1] + 1, S.xi[e][2] + 1
        jN, jC = findfirst(==(iN), obs), findfirst(==(iC), obs)
        isapprox(c3.x0nom[jN], nom.N0; rtol=1e-6) && isapprox(c3.x0nom[jC], nom.C0; rtol=1e-6) && return (e=e, iN=iN, iC=iC)
    end
    error("no $(treated ? "treated" : "untreated") Stage 3a coculture series with ratio $ratio")
end
# Its fitted starting counts (the carried Stage 3a model's u0: estimated where the fit estimates them, nominal otherwise)
function fitted_start(ratio::AbstractString, treated::Bool)
    c3 = zgroup("s3").cands[1]; s = coculture_env(ratio, treated)
    (N0=c3.u0[s.iN], C0=c3.u0[s.iC])
end
# Starting counts of every control run (control_config.toml [initial] source, user 2026-11-12): "fitted_treated"
# (default) starts each ratio from the Stage 3a fitted start of its treated (Ce1) series, so the controller and
# the fit simulate the same wells; "nominal" = the nominal seeding (before 2026-11-12).
initial_source() = String(get(CFG["initial"], "source", "nominal"))
start_counts(ratio::AbstractString) = initial_source() == "fitted_treated" ? fitted_start(ratio, true) : nominal_seeding(ratio)
# measured day-14 counts (N, C) of a Stage 3a coculture series and the fit's residual RMSE over that series' points
function coculture_day14(ratio::AbstractString, treated::Bool)
    g3 = zgroup("s3"); c3 = g3.cands[1]; s = coculture_env(ratio, treated)
    sol = solve(ODEProblem(c3.model, c3.u0, (0.0, 14.0), c3.params), Tsit5(); reltol=1e-10, abstol=1e-10)
    dN = only(d for d in g3.datasets if d.state_index == s.iN); dC = only(d for d in g3.datasets if d.state_index == s.iC)
    at14(d) = d.y[findlast(==(14.0), d.x)]
    rmse(d, i) = sqrt(sum(abs2, [sol(t)[i] for t in d.x] .- d.y) / length(d.x))
    (N=at14(dN), C=at14(dC), rmseN=rmse(dN, s.iN), rmseC=rmse(dC, s.iC), fitN=sol(14.0)[s.iN], fitC=sol(14.0)[s.iC])
end
# Plant against data (user, 2026-11-12): each scenario's plant from each condition's own fitted start, untreated
# (u = 0, against Ce0) and at a constant 1 µM (against Ce1), day-14 N and C next to the measured counts
function plant_vs_data(; Δ=Float64(CFG["schedule"]["decision_interval_days"]))
    rows = NamedTuple[]
    for scen in (:A, :B)
        P = scenario_plant(scen)
        for ratio in ("25:75", "50:50", "75:25"), (treated, u) in ((false, 0.0), (true, 1.0))
            s = fitted_start(ratio, treated); x = x0_for(P, s.N0, s.C0)
            for _ in 1:round(Int, 14.0 / Δ); x = F(x, u, P, Δ); end
            d = coculture_day14(ratio, treated)
            push!(rows, (scenario=scen, ratio=ratio, condition=treated ? "Ce1 (1 µM)" : "Ce0 (untreated)", N0=s.N0, C0=s.C0,
                         N=x[P.iN], C=x[P.iC], N_data=d.N, C_data=d.C, rmseN=d.rmseN, rmseC=d.rmseC))
        end
    end
    rows
end

# Scenario roles and the caveats every control report states: control_text.jl (shared with control_sweep_report.jl)
include(joinpath(@__DIR__, "control_text.jl"))
function b_whatif_note(rows=plant_vs_data())
    tot(sc) = (r = only(x for x in rows if x.scenario == sc && x.ratio == "50:50" && startswith(x.condition, "Ce1")); (r.N + r.C, r.N_data + r.C_data))
    (b, d), (a, _) = tot(:B), tot(:A)
    b_whatif_text(b, d, a)
end

# initial state: the given counts, drug and toxicity states at 0
function x0_for(P::Plant, N0, C0)
    x = zeros(P.nx); x[P.iN] = N0; x[P.iC] = C0
    P.fp === nothing || (x[P.fp.iW] = 1.0)
    L = P.sys.lines
    L[2].action == :resource && error("resource states need initialisation")
    x
end
