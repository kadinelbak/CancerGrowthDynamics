# Model selection rules for the zoo (used by run_zoo_fits.jl and zoo_rebuild_pages.jl).
#
# Per fit: among candidates within RULE_BIC of the lowest BIC, the fewest free parameters wins
#   (lower BIC breaks ties).
# Stage 2 (cisplatin): both lines must share one mechanism structure (kernel + effect) with their
#   own rates. A2780cis is derived from A2780 and cisplatin kills both through DNA damage then
#   delayed death; resistance changes rates (uptake, glutathione, repair, damage tolerance), not
#   the mechanism. The rule is applied to the joint ranking: BIC = BIC_A2780 + BIC_A2780cis (each line
#   is fitted separately with its own error variance), k = k_N + k_C.
# Stage 1b (low resource): a shared structure is tested, not assumed. A = each line's own pick;
#   B = the joint pick. B is used only if its joint BIC is within RULE_BIC of A's summed BIC and its
#   joint k is no larger than A's summed k.

function rule_pick(bics, ks)
    m = minimum(b for b in bics if isfinite(b); init=Inf)
    ok = [j for j in eachindex(bics) if isfinite(bics[j]) && bics[j] <= m + RULE_BIC]
    sort(ok; by=j -> (ks[j], bics[j]))[1]
end
function order_rule(rows)
    s = rule_pick([r.bic for r in rows], [length(r.params) for r in rows])
    vcat(rows[s:s], [r for (j, r) in enumerate(rows) if j != s])
end
# put the named candidate first, the rest by BIC
function order_by_name(rows, name)
    s = findfirst(r -> r.name == name, rows)
    rest = sort([r for (j, r) in enumerate(rows) if j != s]; by=r -> r.bic)
    vcat(rows[s:s], rest)
end

# joint ranking of the structures both lines were fitted with
function joint_ranking(rowsN, rowsC)
    dC = Dict(r.name => r for r in rowsC)
    out = [(name=r.name, label=r.label, bicN=r.bic, bicC=dC[r.name].bic, bic=r.bic + dC[r.name].bic,
            kN=length(r.params), kC=length(dC[r.name].params), k=length(r.params) + length(dC[r.name].params))
           for r in rowsN if haskey(dC, r.name)]
    sort!(out; by=x -> x.bic)
end
joint_pick(jr) = jr[rule_pick([x.bic for x in jr], [x.k for x in jr])]

# Stage 2: returns (pick name, info); rows of both lines get the shared structure first
function select_shared(rowsN, rowsC)
    jr = joint_ranking(rowsN, rowsC); p = joint_pick(jr)
    (name=p.name, info=(rule="shared", pick=p, ranking=jr))
end

# Stage 1b: returns (names per line, info with both options)
function select_ab(rowsN, rowsC)
    aN = rowsN[rule_pick([r.bic for r in rowsN], [length(r.params) for r in rowsN])]
    aC = rowsC[rule_pick([r.bic for r in rowsC], [length(r.params) for r in rowsC])]
    A = (nameN=aN.name, nameC=aC.name, labelN=aN.label, labelC=aC.label, bic=aN.bic + aC.bic, k=length(aN.params) + length(aC.params))
    jr = joint_ranking(rowsN, rowsC); B = joint_pick(jr)
    useB = B.bic <= A.bic + RULE_BIC && B.k <= A.k
    (names=useB ? (B.name, B.name) : (aN.name, aC.name), info=(rule="ab", A=A, B=B, useB=useB, dbic=B.bic - A.bic, ranking=jr))
end

# Replays the rules on every bootstrap replicate, pairing replicate b of the two lines.
# gN, gC: the two groups; rN, rC: their bootstrap results. Returns (picksN, picksC) as candidate indices.
function replay_joint(gN, gC, rN, rC, rule)
    B = min(size(rN.bic, 1), size(rC.bic, 1))
    B == 0 && return (Int[], Int[])
    namesN = [c.name for c in gN.cands]; namesC = [c.name for c in gC.cands]
    kN = [length(c.params) for c in gN.cands]; kC = [length(c.params) for c in gC.cands]
    common = [n for n in namesN if n in namesC]
    iN = Dict(n => findfirst(==(n), namesN) for n in common); iC = Dict(n => findfirst(==(n), namesC) for n in common)
    pN = Int[]; pC = Int[]
    for b in 1:B
        jb = [rN.bic[b, iN[n]] + rC.bic[b, iC[n]] for n in common]; jk = [kN[iN[n]] + kC[iC[n]] for n in common]
        sj = rule_pick(jb, jk)
        if rule == "shared"
            push!(pN, iN[common[sj]]); push!(pC, iC[common[sj]])
        else
            aN = rule_pick(rN.bic[b, :], kN); aC = rule_pick(rC.bic[b, :], kC)
            useB = jb[sj] <= rN.bic[b, aN] + rC.bic[b, aC] + RULE_BIC && jk[sj] <= kN[aN] + kC[aC]
            push!(pN, useB ? iN[common[sj]] : aN); push!(pC, useB ? iC[common[sj]] : aC)
        end
    end
    pN, pC
end

# replay the paired-line rules on the bootstrap results of the per-line stages that have one
# (stages.toml `joint_rule`; the default matches it for older checkpoints)
function apply_joint_replay!(groups, results; rules=(("s1b", "ab"), ("s2", "shared")))
    for (kind, rule) in rules
        gN = only(g for g in groups if g.id == "$kind-a2780"); gC = only(g for g in groups if g.id == "$kind-a2780cis")
        pN, pC = replay_joint(gN, gC, results[gN.id].res, results[gC.id].res, rule)
        isempty(pN) && continue
        results[gN.id] = merge(results[gN.id], (res=merge(results[gN.id].res, (pick=pN,)),))
        results[gC.id] = merge(results[gC.id], (res=merge(results[gC.id].res, (pick=pC,)),))
    end
    results
end

# the selection info for the report when it was not stored with the fits (older checkpoints)
function joint_info_from_groups(groups)
    c(id) = only(g for g in groups if g.id == id).cands
    s2 = select_shared(c("s2-a2780"), c("s2-a2780cis")); s1b = select_ab(c("s1b-a2780"), c("s1b-a2780cis"))
    for (id, n) in (("s2-a2780", s2.name), ("s2-a2780cis", s2.name), ("s1b-a2780", s1b.names[1]), ("s1b-a2780cis", s1b.names[2]))
        c(id)[1].name == n || @warn "the rules now pick $n for $id, but the fits carried $(c(id)[1].name) forward; rerun with ZOO_RESUME=1"
    end
    Dict("s2" => s2.info, "s1b" => s1b.info)
end
