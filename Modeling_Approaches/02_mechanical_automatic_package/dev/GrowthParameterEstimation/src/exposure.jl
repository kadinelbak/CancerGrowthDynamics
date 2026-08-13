module Exposure

export AbstractExposure, ConstantExposure, PulseExposure, SteppedExposure, DecayingExposure,
       build_exposure, evaluate_exposure

abstract type AbstractExposure end

struct ConstantExposure <: AbstractExposure
    value::Float64
end

(structure::ConstantExposure)(t::Real) = structure.value

struct PulseExposure <: AbstractExposure
    amplitude::Float64
    start_time::Float64
    end_time::Float64
end

function (structure::PulseExposure)(t::Real)
    return (structure.start_time <= t <= structure.end_time) ? structure.amplitude : 0.0
end

struct SteppedExposure <: AbstractExposure
    change_times::Vector{Float64}
    values::Vector{Float64}
end

function (structure::SteppedExposure)(t::Real)
    idx = searchsortedlast(structure.change_times, Float64(t))
    idx = clamp(idx, 1, length(structure.values))
    return structure.values[idx]
end

struct DecayingExposure <: AbstractExposure
    c0::Float64
    decay_rate::Float64
    t0::Float64
end

function (structure::DecayingExposure)(t::Real)
    if t < structure.t0
        return 0.0
    end
    return structure.c0 * exp(-structure.decay_rate * (Float64(t) - structure.t0))
end

function build_exposure(kind::Symbol; kwargs...)
    if kind == :constant
        return ConstantExposure(Float64(get(kwargs, :value, 0.0)))
    elseif kind == :pulse
        return PulseExposure(
            Float64(get(kwargs, :amplitude, 1.0)),
            Float64(get(kwargs, :start_time, 0.0)),
            Float64(get(kwargs, :end_time, 1.0)),
        )
    elseif kind == :stepped
        return SteppedExposure(
            Float64.(get(kwargs, :change_times, [0.0])),
            Float64.(get(kwargs, :values, [0.0])),
        )
    elseif kind == :decay
        return DecayingExposure(
            Float64(get(kwargs, :c0, 1.0)),
            Float64(get(kwargs, :decay_rate, 0.1)),
            Float64(get(kwargs, :t0, 0.0)),
        )
    else
        error("Unknown exposure kind: $kind")
    end
end

function evaluate_exposure(exposure::AbstractExposure, times::AbstractVector{<:Real})
    return [exposure(t) for t in times]
end

end # module Exposure
