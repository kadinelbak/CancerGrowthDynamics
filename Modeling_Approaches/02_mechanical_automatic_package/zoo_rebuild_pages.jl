# Redraws the model zoo report and the site bundle from the checkpoint written by
# run_zoo_fits.jl (outputs/zoo_fits/zoo_state.jls), without refitting anything.
#   julia --project=. zoo_rebuild_pages.jl

include(joinpath(@__DIR__, "zoo_common.jl"))
include(joinpath(@__DIR__, "zoo_models.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization

const ZOO_DIR  = joinpath(@__DIR__, "outputs", "zoo_fits")
const ZOO_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits")
const BLUE, ORANGE, GREEN, PURPLE = "#0b6e99", "#d95f02", "#3b8b4f", "#8b5cf6"
state = deserialize(joinpath(ZOO_DIR, "zoo_state.jls"))
const BOOT_B = state.B
include(joinpath(@__DIR__, "zoo_selection.jl"))
JOINT_INFO = hasproperty(state, :joint) ? state.joint : joint_info_from_groups(state.groups)
apply_joint_replay!(state.groups, state.results)
include(joinpath(@__DIR__, "export_zoo_bundle.jl"))
include(joinpath(@__DIR__, "zoo_report.jl"))
write_zoo_report(state.groups, state.results)
export_zoo_bundle(state.groups, state.results)
