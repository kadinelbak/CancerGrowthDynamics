# Scenario roles and the caveats every control report states (user, 2026-11-12). Needs CFG (control_config.toml).
# Included by control_model.jl and by control_sweep_report.jl, which does not load the plants.
# Dose levels the MPC may choose (user, 2026-11-22): [dose] mpc_levels, the measured doses by default; CONTROL_MPC_LEVELS
# ("0,0.62,1,1.47,2,3") overrides it for a side run. The fixed comparators (constant 1 and 1.47 µM, …) never change.
const MPC_LEVELS = haskey(ENV, "CONTROL_MPC_LEVELS") ? parse.(Float64, split(ENV["CONTROL_MPC_LEVELS"], ",")) :
                   Float64.(get(CFG["dose"], "mpc_levels", CFG["dose"]["levels_data"]))
const MPC_MEASURED = sort(MPC_LEVELS) == sort(Float64.(CFG["dose"]["levels_data"]))
const MPC_CAP = maximum(MPC_LEVELS)
const MPC_D = MPC_MEASURED ? "MPC, measured doses only" : "MPC, doses up to $(MPC_CAP) µM"                  # 14-day, dose only
const MPC_F = MPC_MEASURED ? "MPC, dose + fasting (measured levels)" : "MPC, dose + fasting (doses up to $(MPC_CAP) µM)"   # long run
# Extra MPC policies with a higher dose cap (user, 2026-11-22: "one MPC that goes up to three micromolar and one that
# does not"): [dose] mpc_high_levels (empty by default = no extra policy; CONTROL_MPC_HIGH_LEVELS overrides). They run
# next to the measured-dose MPC with the same objective, so both appear side by side.
const MPC_HIGH_LEVELS = haskey(ENV, "CONTROL_MPC_HIGH_LEVELS") ? parse.(Float64, split(ENV["CONTROL_MPC_HIGH_LEVELS"], ",")) :
                        Float64.(get(CFG["dose"], "mpc_high_levels", Float64[]))
const MPC_HIGH_CAP = isempty(MPC_HIGH_LEVELS) ? NaN : maximum(MPC_HIGH_LEVELS)
const MPC_HIGH_STR = replace(string(MPC_HIGH_CAP), r"\.0$" => "")
const MPC_D_HIGH = "MPC drug-only ≤$(MPC_HIGH_STR) µM"
const MPC_F_HIGH = "MPC drug+fasting ≤$(MPC_HIGH_STR) µM"
const HEADLINE_SCEN = String(get(get(CFG, "scenarios", Dict{String,Any}()), "headline", "A"))
scenario_role(scen) = string(scen) == HEADLINE_SCEN ? "main result" :
    string(scen) == "B" ? "what-if: monoculture drug response in coculture" : "sensitivity case"
# b: scenario B's day-14 total at a constant 1 µM (50:50, fitted Ce1 start); d: measured; a: scenario A's
b_whatif_text(b, d, a) = "Scenario B is a what-if: the monoculture drug response applied in coculture. The treated cocultures contradict it: " *
    "at a constant 1 µM it predicts a day-14 total of $(round(Int, b)) cells at 50:50 against $(round(Int, d)) measured " *
    "(scenario A, the Stage 3a coculture fit: $(round(Int, a)))."
const KILL_LIMIT_NOTE = "The coculture kill factors (φ) were measured at 1 µM only, the one dose the cocultures received, and are applied " *
    "unchanged at every dose, so scenario A at 0.62 and 1.47 µM extrapolates them; scenario B is the bound for a resistant kill as in monoculture."
const EXTRAP_NOTE = "Everything past day 14 is model extrapolation: the experiments end at day 14."
const LOWRES_NOTE = "In low resource (0.5% FBS) the resistant line dies without any drug, in the data and in the plant. The treated " *
    "low-resource plant rests on a single series, the 50:50 combination culture (Stage 4)."
