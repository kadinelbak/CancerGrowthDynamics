module MechanicalAutomaticModeling

include("io_utils.jl")
include("model_registry.jl")
include("fit_workflows.jl")
include("analysis_workflows.jl")
include("staged_a2780_workflow.jl")
include("sample_aware_report.jl")
include("reduced_stage_report.jl")
include("adaptive_simulation_engine.jl")
include("a2780_adaptive_adapter.jl")
include("hybrid_refit_service.jl")

using .IOUtils
using .ModelRegistry
using .FitWorkflows
using .AnalysisWorkflows
using .StagedA2780Workflow
using .SampleAwareReport
using .ReducedStageReport
using .AdaptiveSimulationEngine
using .A2780AdaptiveAdapter
using .HybridRefitService

export IOUtils, ModelRegistry, FitWorkflows, AnalysisWorkflows, StagedA2780Workflow, SampleAwareReport, ReducedStageReport, AdaptiveSimulationEngine, A2780AdaptiveAdapter, HybridRefitService

end
