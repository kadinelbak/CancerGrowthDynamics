# Which saved fit set the current outputs/zoo_fits holds (user, 2026-11-23: "the model used right now for the ongoing
# optimal control should be saved as a given set and new reports should indicate what original fits they use if they
# only redo the optimal control"). A fit set is a frozen copy in outputs/fit_sets/<name>/ with a manifest.toml that
# records its fit run, code commit and the checksum of zoo_state.jls; the live zoo_state.jls is matched against it.
# Standard library only, so the control reports and build_final_report.jl can both include it.

using TOML, SHA

const FIT_SETS_DIR = joinpath(@__DIR__, "outputs", "fit_sets")

function fit_set_info(state=joinpath(@__DIR__, "outputs", "zoo_fits", "zoo_state.jls"))
    isfile(state) || return Dict{String,Any}("name" => "none", "matched" => false)
    h = bytes2hex(open(sha256, state))
    if isdir(FIT_SETS_DIR)
        for d in sort(readdir(FIT_SETS_DIR))
            mf = joinpath(FIT_SETS_DIR, d, "manifest.toml")
            isfile(mf) || continue
            m = TOML.parsefile(mf)
            get(m, "state_sha256", "") == h && return merge(Dict{String,Any}(m), Dict{String,Any}("matched" => true, "state_sha256" => h))
        end
    end
    Dict{String,Any}("name" => "unsaved", "matched" => false, "state_sha256" => h,
                     "fit_run" => "the fits in outputs/zoo_fits (not saved as a fit set)")
end

# one line for the reports: the set's name, its fit run and the code commit of the fits
function fits_used_text(fs=fit_set_info())
    get(fs, "matched", false) || return "Fits used: $(get(fs, "fit_run", "unknown")) (zoo_state.jls sha256 $(first(get(fs, "state_sha256", "?"), 12)))."
    "Fits used: saved fit set $(fs["name"]), fit run $(fs["fit_run"]), fit code commit $(first(string(get(fs, "fit_commit", "")), 7)), " *
    "outputs committed in $(get(fs, "outputs_commit", "?")), first Final Report $(get(fs, "first_final_report", "?"))."
end
