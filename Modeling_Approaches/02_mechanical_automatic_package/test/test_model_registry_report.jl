using GrowthParameterEstimation

@testset "Model registry report contract" begin
    package_root = normpath(joinpath(@__DIR__, ".."))
    report_path = joinpath(package_root, "outputs", "reports", "a2780_model_registry.html")
    @test isfile(report_path)

    report = read(report_path, String)
    @test occursin("A2780 Model Registry", report)
    @test occursin("id=\"selectionTray\"", report)
    @test occursin("id=\"selectionPanel\"", report)
    @test occursin("classList.toggle('hidden', selected.size === 0)", report)
    @test occursin("id=\"groupFilter\"", report)
    @test occursin("id=\"statusFilter\"", report)
    @test occursin("joint_ic_effect_hill_ramp_onset", report)
    @test occursin("strobl_birth_growth_cost_turnover", report)
    @test occursin("load_plus_tolerant_growth_context", report)
    @test occursin("transit_chain_erlang", report)
    @test occursin("Legacy local registry", report)
    @test occursin("Select visible", report)
    @test occursin("Selected equations", report)
    @test occursin("id=\"tagButtons\"", report)
    @test occursin("function modelTags(model)", report)
    @test occursin("matchesAnyTag", report)
    @test occursin("model.stages.includes(stage)", report)
    @test occursin("Simple shortlist", report)
    @test occursin("const simpleCandidates", report)
    @test occursin("Why these models were selected", report)
    @test occursin("Catalog scope", report)
    @test occursin("<summary><span class=\"icon\" aria-hidden=\"true\">&#9881;</span>Advanced</summary>", report)
    @test !occursin("Browse the equations tested", report)
    @test !occursin("Start with Workflow models", report)
    @test occursin("States, parameters, and interpretation", report)
    @test occursin("data-group-count", report)

    registry_ids = [match.captures[1] for match in eachmatch(r"M\(\"[^\"]+\",\"([^\"]+)\"", report)]
    @test length(registry_ids) == 66
    @test allunique(registry_ids)

    shortlist_match = match(r"const simpleCandidates = \[([\s\S]*?)\];", report)
    @test shortlist_match !== nothing
    shortlist_ids = [candidate.captures[1] for candidate in eachmatch(r"'([^']+)'", shortlist_match.captures[1])]
    @test length(shortlist_ids) == 14
    @test allunique(shortlist_ids)
    @test all(id -> id in registry_ids, shortlist_ids)

    for model_name in GrowthParameterEstimation.list_models()
        @test model_name in registry_ids
    end
    for spec in MechanicalAutomaticModeling.ModelRegistry.local_model_specs()
        @test spec.name in registry_ids
    end
    for hypothesis in MechanicalAutomaticModeling.FitWorkflows.LINKED_TREATMENT_HYPOTHESES
        @test hypothesis in registry_ids
    end
    for variant in MechanicalAutomaticModeling.FitWorkflows.STROBL_MODEL_VARIANTS
        @test variant.name in registry_ids
    end
end
