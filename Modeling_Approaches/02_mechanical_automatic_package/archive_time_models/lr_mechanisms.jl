# Candidate mechanisms for the low-resource (0.5% FBS) monoculture effect. Include after
# run_new_data_staged_fits.jl and bootstrap_tools.jl. r_X and K_X are held at their Stage 1
# logistic values; the observed population is always state 1 (hidden states: resource R,
# by-product W, quiescent Q, damaged D1/D2). Only the first three models are on the
# Candidate models page; the rest are new here.

lr_sig(z) = 1 / (1 + exp(-z))

function lr_mechanism_specs(X::String, r::Float64, K::Float64, y0::Float64)
    G = "r_{$X}\\,$X\\left(1-\\frac{$X}{K_{$X}}\\right)"
    lhs = "\\frac{d$X}{dt}="
    Z(p) = zero(p[1])
    grow(x) = r*x*(1 - x/K)
    [
     (name="lr_subtractive_loss", family="Constant loss", label="Constant low-resource loss ★",
      pnames=["d_$(X)^LR"], bounds=[(0.0, 5.0)], p0s=[[0.3], [1.0]], u0=[y0],
      tex=lhs * G * "-d^{LR}_{$X}\\,u_R\\,$X",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = grow(x) - p[1]*U_R*x)),
     (name="lr_rho_modifiers", family="Rescaled growth", label="Growth-rate and K multipliers",
      pnames=["ρ_r", "ρ_K"], bounds=[(0.0, 20.0), (1e-3, 3.0)], p0s=[[1.0, 0.8], [2.0, 0.2]], u0=[y0],
      tex=lhs * "\\rho_r\\,r_{$X}\\,$X\\left(1-\\frac{$X}{\\rho_K K_{$X}}\\right)",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = p[1]*r*x*(1 - x/(p[2]*K)))),
     (name="lr_separate_r", family="Rescaled growth", label="Separately fitted growth rate",
      pnames=["r_$(X)^LR"], bounds=[(1e-4, 5.0)], p0s=[[0.5], [0.05]], u0=[y0],
      tex=lhs * "r^{LR}_{$X}\\,$X\\left(1-\\frac{$X}{K_{$X}}\\right)",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = p[1]*x*(1 - x/K))),
     (name="delayed_loss_exp", family="Delayed loss", label="Loss that ramps up (exponential onset)",
      pnames=["d", "λ"], bounds=[(0.0, 10.0), (1e-3, 10.0)], p0s=[[0.5, 0.3], [1.5, 0.1], [0.2, 2.0]], u0=[y0],
      tex=lhs * G * "-d\\left(1-e^{-\\lambda t}\\right)$X",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = grow(x) - p[1]*(1 - exp(-p[2]*t))*x)),
     (name="delayed_loss_sigmoid", family="Delayed loss", label="Loss that switches on at day τ",
      pnames=["d", "τ"], bounds=[(0.0, 10.0), (0.0, 20.0)], p0s=[[1.0, 5.0], [0.5, 8.0], [2.0, 3.0], [1.5, 4.8], [0.3, 6.0]], u0=[y0],
      tex=lhs * G * "-\\frac{d}{1+e^{-(t-\\tau)}}\\,$X",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = grow(x) - p[1]*lr_sig(t - p[2])*x)),
     (name="delayed_loss_sigmoid_w", family="Delayed loss", label="Loss that switches on at day τ over width w",
      pnames=["d", "τ", "w"], bounds=[(0.0, 10.0), (0.0, 20.0), (0.05, 10.0)], p0s=[[1.0, 5.0, 1.0], [0.8, 7.0, 2.0], [2.0, 4.0, 0.5], [1.5, 4.8, 0.2], [3.0, 5.0, 0.1], [1.0, 3.8, 0.5], [0.3, 6.0, 1.0]], u0=[y0],
      tex=lhs * G * "-\\frac{d}{1+e^{-(t-\\tau)/w}}\\,$X",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = grow(x) - p[1]*lr_sig((t - p[2])/p[3])*x)),
     (name="decaying_K", family="Shrinking capacity", label="Carrying capacity decays over time",
      pnames=["δ"], bounds=[(0.0, 1.0)], p0s=[[0.1], [0.3], [0.03]], u0=[y0],
      tex=lhs * "r_{$X}\\,$X\\left(1-\\frac{$X}{K_{$X}\\,e^{-\\delta t}}\\right)",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = r*x*(1 - x/max(K*exp(-p[1]*t), 1e-6)))),
     (name="decaying_r_loss", family="Fading growth", label="Growth rate fades, constant loss",
      pnames=["δ", "d"], bounds=[(0.0, 5.0), (0.0, 5.0)], p0s=[[0.2, 0.2], [0.5, 0.5], [0.1, 1.0]], u0=[y0],
      tex=lhs * "r_{$X}e^{-\\delta t}\\,$X\\left(1-\\frac{$X}{K_{$X}}\\right)-d\\,$X",
      model=(du, u, p, t) -> (x = max(u[1], Z(p)); du[1] = exp(-p[1]*t)*grow(x) - p[2]*x)),
     (name="resource_depletion", family="Resource depletion", label="Cells use up the resource, then starve",
      pnames=["c", "d"], bounds=[(0.0, 100.0), (0.0, 10.0)], p0s=[[1.0, 1.0], [5.0, 0.5], [0.3, 3.0]], u0=[y0, 1.0],
      tex=lhs * "r_{$X}R\\,$X\\left(1-\\frac{$X}{K_{$X}}\\right)-d\\,(1-R)\\,$X,\\qquad \\frac{dR}{dt}=-c\\,R\\,\\frac{$X}{1000},\\ R(0)=1",
      model=(du, u, p, t) -> begin
          x = max(u[1], Z(p)); R = max(u[2], Z(p))
          du[1] = R*grow(x) - p[2]*(1 - R)*x
          du[2] = -p[1]*R*x/1000
      end),
     (name="toxin", family="By-product build-up", label="Toxic by-product accumulates",
      pnames=["γ"], bounds=[(0.0, 10.0)], p0s=[[0.1], [0.5], [0.02]], u0=[y0, 0.0],
      tex=lhs * G * "-\\gamma\\,W\\,$X,\\qquad \\frac{dW}{dt}=\\frac{$X}{1000},\\ W(0)=0",
      model=(du, u, p, t) -> begin
          x = max(u[1], Z(p))
          du[1] = grow(x) - p[1]*u[2]*x
          du[2] = x/1000
      end),
     (name="toxin_clearance", family="By-product build-up", label="Toxic by-product accumulates and clears",
      pnames=["γ", "κ"], bounds=[(0.0, 10.0), (0.0, 10.0)], p0s=[[0.1, 0.1], [0.5, 0.5], [0.05, 0.01]], u0=[y0, 0.0],
      tex=lhs * G * "-\\gamma\\,W\\,$X,\\qquad \\frac{dW}{dt}=\\frac{$X}{1000}-\\kappa W,\\ W(0)=0",
      model=(du, u, p, t) -> begin
          x = max(u[1], Z(p))
          du[1] = grow(x) - p[1]*u[2]*x
          du[2] = x/1000 - p[2]*u[2]
      end),
     (name="quiescence_death", family="Switch to non-dividing state", label="Cells become quiescent, then die",
      pnames=["s", "δ"], bounds=[(0.0, 10.0), (0.0, 10.0)], p0s=[[0.5, 0.3], [1.0, 0.5], [0.2, 1.0]], u0=[y0, 0.0],
      tex=lhs * "r_{$X}P\\left(1-\\frac{$X}{K_{$X}}\\right)-\\delta\\,Q,\\qquad \\frac{dQ}{dt}=s\\,P-\\delta\\,Q,\\quad P=$X-Q",
      model=(du, u, p, t) -> begin
          T = max(u[1], Z(p)); Q = max(u[2], Z(p)); P = max(T - Q, Z(p))
          du[1] = r*P*(1 - T/K) - p[2]*Q
          du[2] = p[1]*P - p[2]*Q
      end),
     (name="transit_damage", family="Switch to non-dividing state", label="Cells are damaged, then die through two stages",
      pnames=["k", "τ_D"], bounds=[(0.0, 10.0), (1e-3, 10.0)], p0s=[[0.8, 0.5], [1.5, 1.0], [0.5, 0.3]], u0=[y0, 0.0, 0.0],
      tex=lhs * "r_{$X}P\\left(1-\\frac{$X}{K_{$X}}\\right)-\\tau_D D_2,\\quad \\frac{dD_1}{dt}=k\\,P-\\tau_D D_1,\\quad \\frac{dD_2}{dt}=\\tau_D(D_1-D_2),\\quad P=$X-D_1-D_2",
      model=(du, u, p, t) -> begin
          T = max(u[1], Z(p)); D1 = max(u[2], Z(p)); D2 = max(u[3], Z(p)); P = max(T - D1 - D2, Z(p))
          du[1] = r*P*(1 - T/K) - p[2]*D2
          du[2] = p[1]*P - p[2]*D1
          du[3] = p[2]*(D1 - D2)
      end),
    ]
end

# Fit every mechanism from several starts; keep the lowest BIC. Returns rows sorted by BIC.
function fit_lr_mechanisms(line)
    s1 = stage1[line]; s = stage1b[line]
    datasets = [(x=s.x[.!s.bad], y=s.y[.!s.bad], state_index=1)]
    rows = NamedTuple[]
    for c in lr_mechanism_specs(s1.sym, s1.r, s1.K, s.y[1])
        best = (params=fill(NaN, length(c.pnames)), bic=Inf, sse=Inf)
        # generous time limit: results must not depend on how busy the machine is
        for p0 in c.p0s
            f = bs_fit(c.model, datasets, c.u0, p0, c.bounds; maxtime=900.0)
            f.bic < best.bic && (best = f)
        end
        isfinite(best.bic) && (f = bs_fit(c.model, datasets, c.u0, best.params, c.bounds; maxtime=900.0); f.bic < best.bic && (best = f))
        isfinite(best.bic) || (@warn "no fit" line c.name; continue)
        se = try param_se(c.model, datasets, c.u0, best.params) catch; fill(NaN, length(best.params)) end
        sol = solve(ODEProblem(c.model, c.u0, (0.0, 200.0), best.params), Tsit5(); reltol=1e-8, abstol=1e-8)
        x30, x60, x200 = sol(30.0)[1], sol(60.0)[1], sol(200.0)[1]
        push!(rows, (name=c.name, label=c.label, family=c.family, model=c.model, pnames=c.pnames, params=best.params,
                     bic=best.bic, sse=best.sse, optimizer=:multistart, bounds=c.bounds, se=se, u0=c.u0, tex=c.tex,
                     x30=x30, x60=x60, x200=x200, extinct=x200 < 1.0))
    end
    # Stage 1b candidates added in run_new_data_staged_fits.jl that are not in the list above
    # (accumulated-damage models); polished from their fitted values, same long-term check.
    mine = Set(r.name for r in rows)
    G = "r_{$(s1.sym)}\\,$(s1.sym)\\left(1-\\frac{$(s1.sym)}{K_{$(s1.sym)}}\\right)"; lhs = "\\frac{d$(s1.sym)}{dt}="
    extra_tex = Dict(
        "lr_accumulated_loss4" => lhs * G * "-d^{LR}\\,u_R\\,E_4(\\lambda t)\\,$(s1.sym)",
        "lr_loss_plus_accum4" => lhs * G * "-\\left(d^{LR}+d^{acc}E_4(\\lambda t)\\right)u_R\\,$(s1.sym)",
        "lr_loss_plus_accum8" => lhs * G * "-\\left(d^{LR}+d^{acc}E_8(\\lambda t)\\right)u_R\\,$(s1.sym)")
    for c in s.rows
        c.name in mine && continue
        u0 = [s.y[1]]
        f = bs_fit(c.model, datasets, u0, c.params, c.bounds; maxtime=900.0)
        p, bic, sse = f.bic < c.bic ? (f.params, f.bic, f.sse) : (c.params, c.bic, c.sse)
        se = try param_se(c.model, datasets, u0, p) catch; fill(NaN, length(p)) end
        sol = solve(ODEProblem(c.model, u0, (0.0, 200.0), p), Tsit5(); reltol=1e-8, abstol=1e-8)
        push!(rows, (name=c.name, label=c.label, family="Accumulated damage", model=c.model, pnames=c.pnames, params=p,
                     bic=bic, sse=sse, optimizer=:multistart, bounds=c.bounds, se=se, u0=u0,
                     tex=get(extra_tex, c.name, c.label), x30=sol(30.0)[1], x60=sol(60.0)[1], x200=sol(200.0)[1], extinct=sol(200.0)[1] < 1.0))
    end
    sort!(rows; by=r -> r.bic)
end
