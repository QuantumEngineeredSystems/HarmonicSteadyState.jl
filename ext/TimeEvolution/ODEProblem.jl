using OrdinaryDiffEqTsit5: ODEProblem, solve, ODESolution

# Accept both the time-dependent harmonic symbols and the bracket-free symbols
# returned by get_single_solution, preserving the equation's variable order.
function initial_values(eom::HarmonicEquation, values)
    return map(zip(get_variables(eom), declare_variables(eom))) do (variable, bare)
        haskey(values, variable) && return values[variable]
        haskey(values, bare) && return values[bare]
        return throw(ArgumentError("No initial value for $bare; pass it or set `u0`."))
    end
end

"""
    $(TYPEDSIGNATURES)

Creates an ODEProblem object used by OrdinaryDiffEqTsit5.jl from the equations in `eom` to simulate time-evolution within `timespan`.
`fixed_parameters` must be a dictionary mapping parameters+variables to numbers (possible to use a solution index, e.g. solutions[x][y] for branch y of solution x).
The right-hand side is compiled once; each evaluation only computes numerical sweep values.
If `u0` is specified, it is used as an initial condition; otherwise state values are read in
harmonic-variable order from `fixed_parameters`. Both harmonic symbols and the bracket-free
symbols returned by `get_single_solution` are accepted.
"""
function OrdinaryDiffEqTsit5.ODEProblem(
    eom::HarmonicEquation,
    fixed_parameters;
    sweep::AdiabaticSweep=AdiabaticSweep(),
    u0::Vector=[],
    timespan::Tuple,
    perturb_initial=0.0,
    kwargs...,
)
    if !is_rearranged(eom) # check if time-derivatives of the variable are on the right hand side
        eom = rearrange_standard(eom)
    end

    vars = get_variables(eom)
    fixed = filter_duplicate_parameters(sweep, Dict(fixed_parameters))
    # State values initialize the trajectory; they must not freeze its variables.
    for variable in vcat(vars, declare_variables(eom))
        delete!(fixed, variable)
    end
    rhs = real.(substitute_all(Num.(getfield.(eom.equations, :lhs)), fixed))
    parameters = collect(keys(sweep))
    functions = Tuple(sweep[parameter] for parameter in parameters)
    time = only(get_independent_variables(eom))
    _, compiled! = Symbolics.build_function(
        rhs, vars, parameters, time; expression=Val(false)
    )

    function f!(du, u, p, T)
        values = map(functions) do function_at_time
            return function_at_time(T)
        end
        compiled!(du, u, values, T)
        return nothing
    end

    initial = if isempty(u0)
        real.(initial_values(eom, fixed_parameters)) * (1 - perturb_initial)
    else
        u0
    end

    return ODEProblem(f!, initial, timespan; kwargs...)
end

function HarmonicSteadyState.transform_solutions(
    soln::ODESolution, f::String, harm_eq::HarmonicEquation
)
    return transform_solutions(soln.u, f, harm_eq)
end
function HarmonicSteadyState.transform_solutions(
    s::ODESolution, funcs::Vector{String}, harm_eq::HarmonicEquation
)
    return [transform_solutions(s, f, harm_eq) for f in funcs]
end
