module ModelRegistry

using GrowthParameterEstimation

export ensure_model_registry!, register_custom_models!

const _registry_loaded = Ref(false)

function register_custom_models!()
    # Placeholder hook for full custom model zoo integration.
    # Keep this function idempotent so notebooks can safely rerun cells.
    # Add register_model(ModelSpec(...)) calls here as custom models are ported.
    return nothing
end

function ensure_model_registry!()
    _registry_loaded[] && return nothing
    # Built-ins are available from GrowthParameterEstimation by default.
    register_custom_models!()
    _registry_loaded[] = true
    return nothing
end

end
