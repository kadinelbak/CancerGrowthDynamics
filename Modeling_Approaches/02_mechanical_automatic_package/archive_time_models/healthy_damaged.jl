# Healthy / damaged two-population model, shared by run_new_data_staged_fits.jl and
# run_new_data_validation.jl.
#
# Each line X (N = A2780 naive/sensitive, C = A2780cis resistant) has healthy cells H_X that
# divide, and damaged cells that do not divide but still take up space and are still counted.
# A stress (cisplatin or low resource) is a damage "channel" c with rate k_c; damaged cells
# pass through n_c compartments at rate δ_c and die when they leave the last one, or are
# repaired from the first one at rate ρ_c:
#   dH_X/dt    = r_X H_X (1 − (X + α_XY Y)/K) − Σ_c k_c H_X + Σ_c ρ_c D_c,1
#   dD_c,1/dt  = k_c H_X − (ρ_c + δ_c) D_c,1
#   dD_c,j/dt  = δ_c (D_c,j−1 − D_c,j),   j = 2..n_c
#   X          = H_X + Σ D   (the observed count)
# Literature: damaged-cell transit compartments (Simeoni et al. 2004, Cancer Res 64:1094;
# Lobo & Balthasar 2002, AAPS PharmSci 4:E42), repair of cisplatin lesions that is faster in
# resistant A2780 derivatives (Masuda et al. 1988, Cancer Res 48:5713; Parker et al. 1991,
# J Clin Invest 87:772), and reversible proliferating/quiescent switching under nutrient
# limitation (Gyllenberg & Webb 1989, Growth Dev Aging 53:25).
#
# States are stored as [X, D_1.., D_n] per line so the observed count is a single state
# (H = X − ΣD); per environment the lines follow each other: [N, D_N.., C, D_C..].

hd_width(ns) = 1 + sum(ns; init=0)          # ns = compartments per channel for one line

# offsets[e, l] + 1 is the observed-count state of line l in environment e
function hd_offsets(nE, ns)
    w = [hd_width(n) for n in ns]; per = sum(w)
    [(e - 1)*per + sum(w[1:l-1]; init=0) for e in 1:nE, l in eachindex(ns)]
end
hd_obs(nE, ns) = vec(permutedims(hd_offsets(nE, ns) .+ 1))   # order: env 1 lines…, env 2 lines…

function hd_u0(nE, ns, X0)                  # X0[e][l] = starting count; damaged start at 0
    off = hd_offsets(nE, ns)
    u = zeros(off[end] + hd_width(ns[end]))
    for e in 1:nE, l in eachindex(ns)
        u[off[e, l] + 1] = X0[e][l]
    end
    u
end

# getp(p) returns (r, K, α, ch, act):
#   r[l], K[l]; α[l][m] = crowding of line l by line m (α[l][l] is ignored, taken as 1);
#   ch[l][c] = (κ, ρ, δ); act(e, l, c, t, g) = activity of channel c in env e at time t, where
#   g = (X + α Y)/K is the crowding fraction (k_c = κ·act).
function make_hd(nE, ns, getp)
    off = hd_offsets(nE, ns); nL = length(ns)
    function f!(du, u, p, t)
        q = getp(p)
        z = zero(eltype(u)) * zero(eltype(p))
        for e in 1:nE, l in 1:nL
            o = off[e, l]
            X = max(u[o+1], z)
            crowd = X
            for m in 1:nL
                m == l || (crowd += q.α[l][m] * max(u[off[e, m]+1], z))
            end
            Dsum = z
            for j in 1:sum(ns[l]; init=0)
                Dsum += max(u[o+1+j], z)
            end
            H = max(X - Dsum, z)
            dX = q.r[l] * H * (1 - crowd / q.K[l])
            i = o + 1
            for (c, n) in enumerate(ns[l])
                κ, ρ, δ = q.ch[l][c]
                k = κ * q.act(e, l, c, t, crowd / q.K[l])
                du[i+1] = k*H - (ρ + δ)*u[i+1]
                for j in 2:n
                    du[i+j] = δ*(u[i+j-1] - u[i+j])
                end
                dX -= δ*u[i+n]
                i += n
            end
            du[o+1] = dX
        end
    end
end
