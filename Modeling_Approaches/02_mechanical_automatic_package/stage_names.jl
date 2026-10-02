# Display names of the fitting stages, as used in the thesis (user, 2026-10-02).
# Internal group ids (s1, s1b, s2, s3, s3b, s4), checkpoints and CSV columns keep the old numbering;
# every page and bundle is passed through stage_names() just before it is written.
#
#   code / checkpoint    1     1b    2     3     3b    4
#   displayed            1     2a    2b    3a    3b    4
#
# "Stage 2", "Stages 1b/3b", "Stages 1, 1b and 2", "Stage 1 / 1b / 2" are all rewritten. The rewrite
# is idempotent: names already in the new form (2a, 2b, 3a) are left as they are.
const STAGE_DISPLAY = Dict("1" => "1", "1b" => "2a", "2" => "2b", "3" => "3a", "3b" => "3b", "4" => "4")
const STAGE_RX = r"\b(Stages?)(\s+)(\d[ab]?(?:\s*(?:,|/|and|or|–)\s*\d[ab]?)*)(?!\w)"

# a group's short stage label from the checkpoint, "1b Low-resource monoculture" → "2a Low-resource monoculture"
stage_label(s::AbstractString) = (m = match(r"^(\d[ab]?)(?!\w)(.*)$"s, s); m === nothing ? s : get(STAGE_DISPLAY, m[1], m[1]) * m[2])

function stage_names(s::AbstractString)
    replace(s, STAGE_RX => function (m)
        c = match(STAGE_RX, m)
        c[1] * c[2] * replace(c[3], r"\d[ab]?" => t -> get(STAGE_DISPLAY, t, t))
    end)
end
