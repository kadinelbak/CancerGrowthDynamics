module Observation

export ObservationSpec, observed_signal, viable_total, sum_states

struct ObservationSpec
    name::String
    map_fn::Function
    scale::Float64
    sigma_add::Float64
end

function observed_signal(spec::ObservationSpec, state::AbstractVector, p, t)
    return spec.scale * spec.map_fn(state, p, t)
end

viable_total(state::AbstractVector, p, t) = state[1]
sum_states(indices::Vector{Int}) = (state, p, t) -> sum(state[i] for i in indices)

end # module Observation
