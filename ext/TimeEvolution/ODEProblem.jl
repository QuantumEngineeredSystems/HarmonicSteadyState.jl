using OrdinaryDiffEqTsit5: ODEProblem, solve, ODESolution, Tsit5
using LinearAlgebra: norm
using HarmonicSteadyState: is_real
using SciMLBase: successful_retcode, remake

# Accept both the time-dependent harmonic symbols and the bracket-free symbols
# returned by get_single_solution, preserving the equation's variable order.
function initial_values(eom::HarmonicEquation, values)
    return [
        haskey(values, variable) ? values[variable] : values[bare] for
        (variable, bare) in zip(get_variables(eom), declare_variables(eom))
    ]
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

"""
$(TYPEDSIGNATURES)

Numerically investigate the stability of a solution `soln` of `eom` within `timespan`.
Perturb each state coordinate in both directions by
`perturb_initial * max(norm(reference), 1)`, including when the reference is zero.
Return `true` when all final states remain within `tol * max(norm(reference), 1)`
of the unperturbed reference. The default `algorithm` is `Tsit5()`.
This is a finite-time numerical check, not a proof of asymptotic stability.

"""
function HarmonicSteadyState.is_stable(
    steady_solution::StateDict,
    eom::HarmonicEquation;
    timespan,
    tol=1e-1,
    perturb_initial=1e-3,
    algorithm=Tsit5(),
)
    reference = initial_values(eom, steady_solution)
    all(is_real, reference) || error("the solution is complex!")
    reference = real.(reference)
    scale = max(norm(reference), 1)
    # Test both signs of every coordinate: one shared direction can lie entirely
    # in a stable subspace of an unstable reference state.
    problem = ODEProblem(eom, steady_solution; u0=reference, timespan)
    for coordinate in eachindex(reference), direction in (-1, 1)
        initial = copy(reference)
        initial[coordinate] += direction * perturb_initial * scale
        solution = solve(remake(problem; u0=initial), algorithm)
        successful_retcode(solution) || return false
        final = solution.u[end]
        all(is_real, final) || error("the solution is complex!")
        norm(final - reference) / scale < tol || return false
    end
    return true
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
