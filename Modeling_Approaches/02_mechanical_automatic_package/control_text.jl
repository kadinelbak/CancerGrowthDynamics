# Scenario roles and the caveats every control report states (user, 2026-11-12). Needs CFG (control_config.toml).
# Included by control_model.jl and by control_sweep_report.jl, which does not load the plants.
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
