module AdaptiveSimulationEngine

using Statistics

export TreatmentWindow,
       TreatmentProtocol,
       SimulationScenario,
       ProtocolResult,
       simulate_protocol,
       simulate_ensemble,
       summarize_outcomes,
       equivalent_untreated_time

struct TreatmentWindow
    start_day::Float64
    stop_day::Float64
    dose::Float64
    refresh_interval::Float64
    label::String

    function TreatmentWindow(start_day, stop_day, dose, refresh_interval; label = "")
        start = Float64(start_day)
        stop = Float64(stop_day)
        amount = Float64(dose)
        refresh = Float64(refresh_interval)
        all(isfinite, (start, stop, amount, refresh)) || throw(ArgumentError("treatment-window values must be finite"))
        start >= 0 || throw(ArgumentError("start_day must be nonnegative"))
        stop > start || throw(ArgumentError("stop_day must be greater than start_day"))
        amount >= 0 || throw(ArgumentError("dose must be nonnegative"))
        refresh > 0 || throw(ArgumentError("refresh_interval must be positive"))
        new(start, stop, amount, refresh, String(label))
    end
end

struct TreatmentProtocol{P}
    windows::Vector{TreatmentWindow}
    monitoring_interval::Float64
    decision_policy::P
    allowable_doses::Vector{Float64}
    horizon::Float64

    function TreatmentProtocol(
        windows;
        monitoring_interval = 1.0,
        decision_policy = nothing,
        allowable_doses = Float64[],
        horizon = nothing,
    )
        ordered = sort!(collect(TreatmentWindow, windows); by = window -> (window.start_day, window.stop_day))
        interval = Float64(monitoring_interval)
        isfinite(interval) && interval > 0 || throw(ArgumentError("monitoring_interval must be finite and positive"))
        doses = sort!(unique!(Float64.(allowable_doses)))
        all(dose -> isfinite(dose) && dose >= 0, doses) || throw(ArgumentError("allowable doses must be finite and nonnegative"))
        active_windows = filter(window -> window.dose > 0, ordered)
        for (left, right) in zip(active_windows, Iterators.drop(active_windows, 1))
            if right.start_day < left.stop_day
                throw(ArgumentError("active treatment windows overlap: $(left.label) and $(right.label)"))
            end
        end
        inferred_horizon = isempty(ordered) ? interval : maximum(window.stop_day for window in ordered)
        final_horizon = horizon === nothing ? inferred_horizon : Float64(horizon)
        isfinite(final_horizon) && final_horizon > 0 || throw(ArgumentError("horizon must be finite and positive"))
        new{typeof(decision_policy)}(ordered, interval, decision_policy, doses, final_horizon)
    end
end

struct SimulationScenario{M,P,S,O,F}
    model::M
    parameters::P
    initial_state::S
    protocol::TreatmentProtocol
    observable::O
    segment_simulator::F
    save_interval::Float64

    function SimulationScenario(
        model,
        parameters,
        initial_state,
        protocol::TreatmentProtocol;
        observable = identity,
        segment_simulator,
        save_interval = 0.1,
    )
        interval = Float64(save_interval)
        isfinite(interval) && interval > 0 || throw(ArgumentError("save_interval must be finite and positive"))
        new{typeof(model),typeof(parameters),typeof(initial_state),typeof(observable),typeof(segment_simulator)}(
            model,
            parameters,
            initial_state,
            protocol,
            observable,
            segment_simulator,
            interval,
        )
    end
end

struct ProtocolResult{S,O}
    times::Vector{Float64}
    states::Vector{S}
    observables::Vector{O}
    commanded_dose::Vector{Float64}
    effective_exposure::Vector{Float64}
    damage_signal::Vector{Float64}
    pharmacodynamic_active::BitVector
    decisions::Vector{NamedTuple}
end

function _window_dose(windows, time)
    for window in windows
        window.start_day <= time < window.stop_day && return window.dose
    end
    return 0.0
end

function _event_times(protocol, save_interval)
    events = collect(0.0:save_interval:protocol.horizon)
    append!(events, 0.0:protocol.monitoring_interval:protocol.horizon)
    for window in protocol.windows
        append!(events, (window.start_day, window.stop_day))
        window.dose > 0 && append!(events, window.start_day:window.refresh_interval:window.stop_day)
    end
    push!(events, protocol.horizon)
    return sort!(unique!(clamp.(events, 0.0, protocol.horizon)))
end

function _policy_dose(protocol, time, state, previous_dose, history)
    protocol.decision_policy === nothing && return _window_dose(protocol.windows, time)
    decision = protocol.decision_policy((
        time = time,
        state = state,
        previous_dose = previous_dose,
        history = history,
        allowable_doses = protocol.allowable_doses,
    ))
    dose = Float64(decision)
    isfinite(dose) && dose >= 0 || throw(ArgumentError("decision policy returned an invalid dose"))
    isempty(protocol.allowable_doses) || any(==(dose), protocol.allowable_doses) ||
        throw(ArgumentError("decision policy returned dose $dose outside allowable_doses"))
    return dose
end

"""
Simulate a protocol by calling `segment_simulator` between every treatment,
refresh, monitoring, and output boundary. The callback receives
`(model, state, parameters, t0, t1, dose)` and returns either a state or a
named tuple with `state` and optional `effective_exposure`/`damage_signal`.
"""
function simulate_protocol(scenario::SimulationScenario)
    times = _event_times(scenario.protocol, scenario.save_interval)
    state = deepcopy(scenario.initial_state)
    states = typeof(state)[deepcopy(state)]
    observables = [scenario.observable(state)]
    commanded = Float64[0.0]
    exposure = Float64[0.0]
    damage = Float64[0.0]
    decisions = NamedTuple[]
    previous_dose = 0.0
    last_monitor = -Inf

    for (t0, t1) in zip(times, Iterators.drop(times, 1))
        monitoring = isapprox(mod(t0, scenario.protocol.monitoring_interval), 0.0; atol = 1e-8)
        dose = if monitoring && t0 > last_monitor + 1e-8
            last_monitor = t0
            selected = _policy_dose(scenario.protocol, t0, state, previous_dose, observables)
            push!(decisions, (time = t0, previous_dose = previous_dose, selected_dose = selected))
            selected
        elseif scenario.protocol.decision_policy === nothing
            _window_dose(scenario.protocol.windows, t0)
        else
            previous_dose
        end
        segment = scenario.segment_simulator(scenario.model, state, scenario.parameters, t0, t1, dose)
        if segment isa NamedTuple
            haskey(segment, :state) || throw(ArgumentError("segment result must contain state"))
            state = segment.state
            effective = Float64(get(segment, :effective_exposure, dose))
            signal = Float64(get(segment, :damage_signal, effective))
        else
            state = segment
            effective = dose
            signal = dose
        end
        all(isfinite, state) || throw(ArgumentError("segment simulation produced a non-finite state"))
        any(value -> value < -sqrt(eps(Float64)), state) && throw(ArgumentError("segment simulation produced a negative state"))
        state = max.(state, zero(eltype(state)))
        push!(states, deepcopy(state))
        push!(observables, scenario.observable(state))
        push!(commanded, dose)
        push!(exposure, effective)
        push!(damage, signal)
        previous_dose = dose
    end
    active_threshold = max(maximum(exposure; init = 0.0) * 0.01, eps(Float64))
    return ProtocolResult(times, states, observables, commanded, exposure, damage, BitVector(exposure .> active_threshold), decisions)
end

function simulate_ensemble(scenario::SimulationScenario, parameter_sets)
    return [simulate_protocol(SimulationScenario(
        scenario.model,
        parameters,
        scenario.initial_state,
        scenario.protocol;
        observable = scenario.observable,
        segment_simulator = scenario.segment_simulator,
        save_interval = scenario.save_interval,
    )) for parameters in parameter_sets]
end

_trapz(times, values) = sum((times[index + 1] - times[index]) * (values[index + 1] + values[index]) / 2 for index in 1:(length(times) - 1); init = 0.0)
_step_auc(times, values) = sum((times[index + 1] - times[index]) * values[index + 1] for index in 1:(length(times) - 1); init = 0.0)

function _observable_total(observable)
    if observable isa Number
        return Float64(observable)
    elseif observable isa NamedTuple && :total in keys(observable)
        return Float64(observable.total)
    elseif observable isa AbstractVector && length(observable) == 3 && isapprox(observable[3], observable[1] + observable[2])
        return Float64(observable[3])
    end
    return Float64(sum(observable))
end

function _observable_cis(observable)
    observable isa Number && return Float64(observable)
    observable isa NamedTuple && :cis in keys(observable) && return Float64(observable.cis)
    return length(observable) >= 2 ? Float64(observable[2]) : Float64(observable[end])
end

function summarize_outcomes(results, progression_rule)
    collection = results isa ProtocolResult ? [results] : collect(results)
    isempty(collection) && throw(ArgumentError("at least one simulation result is required"))
    summaries = map(collection) do result
        total = [_observable_total(observable) for observable in result.observables]
        cis = [_observable_cis(observable) for observable in result.observables]
        progressed = findfirst(index -> progression_rule(result.times[index], result.observables[index]), eachindex(result.times))
        off_time = sum(result.times[index + 1] - result.times[index] for index in 1:(length(result.times) - 1) if result.commanded_dose[index + 1] == 0; init = 0.0)
        cycle_count = count(index -> result.commanded_dose[index] <= 0 && result.commanded_dose[index + 1] > 0, 1:(length(result.commanded_dose) - 1))
        threshold_overshoot = maximum((progression_rule(time, observable) ? _observable_total(observable) : 0.0) for (time, observable) in zip(result.times, result.observables); init = 0.0)
        (
            time_to_progression = progressed === nothing ? Inf : result.times[progressed],
            controlled = progressed === nothing,
            total_auc = _trapz(result.times, total),
            cis_auc = _trapz(result.times, cis),
            final_cis = cis[end],
            final_total = total[end],
            cumulative_dose = _step_auc(result.times, result.commanded_dose),
            exposure_auc = _trapz(result.times, result.effective_exposure),
            time_off_treatment = off_time,
            cycle_count = cycle_count,
            decision_count = length(result.decisions),
            threshold_overshoot = threshold_overshoot,
        )
    end
    return (
        members = summaries,
        probability_of_control = mean(summary.controlled for summary in summaries),
        median_final_total = median(summary.final_total for summary in summaries),
        median_final_cis = median(summary.final_cis for summary in summaries),
    )
end

function equivalent_untreated_time(state, untreated_reference)
    isempty(untreated_reference) && return nothing
    target = sum(state)
    index = argmin(abs(sum(entry.state) - target) for entry in untreated_reference)
    entry = untreated_reference[index]
    return (time = entry.time, distance = abs(sum(entry.state) - target), reference_state = entry.state)
end

end
