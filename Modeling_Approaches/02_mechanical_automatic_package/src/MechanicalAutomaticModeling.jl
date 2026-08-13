module MechanicalAutomaticModeling

include("io_utils.jl")
include("model_registry.jl")
include("fit_workflows.jl")
include("analysis_workflows.jl")
include("staged_a2780_workflow.jl")

using .IOUtils
using .ModelRegistry
using .FitWorkflows
using .AnalysisWorkflows
using .StagedA2780Workflow

export IOUtils, ModelRegistry, FitWorkflows, AnalysisWorkflows, StagedA2780Workflow

end
