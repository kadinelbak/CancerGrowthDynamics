module MechanicalAutomaticModeling

include("io_utils.jl")
include("model_registry.jl")
include("fit_workflows.jl")
include("analysis_workflows.jl")

using .IOUtils
using .ModelRegistry
using .FitWorkflows
using .AnalysisWorkflows

export IOUtils, ModelRegistry, FitWorkflows, AnalysisWorkflows

end
