module HarmonicSteadyState

using QuestBase:
    QuestBase,
    HarmonicEquation,
    var_name,
    declare_variables,
    declare_variable,
    hasnan,
    substitute_all,
    get_independent_variables,
    d,
    _remove_brackets,
    source,
    source_type

# default global settings
IM_TOL::Float64 = 1e-6
function set_imaginary_tolerance(x::Float64)
    @eval(IM_TOL::Float64 = $x)
end

using DocStringExtensions: TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using OrderedCollections: OrderedDict
using ProgressMeter: ProgressMeter, Progress
using LinearAlgebra: LinearAlgebra, eigvals
using Random: Random # for setting seed

using Distances: Distances
using BijectiveHilbert: BijectiveHilbert, Simple2D, decode_hilbert!, encode_hilbert
using HomotopyContinuation: HomotopyContinuation
using Symbolics: Symbolics, unwrap, wrap, Num, get_variables
using SymbolicUtils: SymbolicUtils

const HC = HomotopyContinuation
import FunctionWrappers: FunctionWrapper
using RuntimeGeneratedFunctions: RuntimeGeneratedFunction

include("extension_functions.jl")

include("HC_wrapper.jl")
using .HC_wrapper

include("types.jl")
include("utils.jl")
include("Problem.jl")
include("Jacobian.jl")
include("Result.jl")
include("methods.jl")

include("solve_homotopy.jl")
include("sorting.jl")
include("classification.jl")
include("transform_solutions.jl")

include("LinearResponse/LinearResponse.jl")
using .LinearResponse

include("LimitCycles/LimitCycles.jl")
using .LimitCycles

# Equation
export HarmonicEquation # for Meanfield equations

# methods
export WarmUp
export TotalDegree
export Polyhedral

# handle solutions
export get_steady_states
export classify_solutions!
export get_class
export filter_result!
export get_single_solution
export get_solutions
export get_branches
export transform_solutions
export IM_TOL
export set_imaginary_tolerance

# Result
export swept_parameter, swept_parameters
export attractors
export phase_diagram

# Limit cycles
export get_cycle_variables, get_limit_cycles, add_pairs!

# LinearResponse
export eigenvalues, eigenvectors
export get_jacobian_response
export get_linear_response
export get_rotframe_jacobian_response
export get_susceptibility
export get_forward_transmission_response

# plotting
export plot_linear_response
export plot_phase_diagram
export plot_rotframe_jacobian_response
export plot_eigenvalues
export plot_spaghetti

# extension functions
export AdiabaticSweep
export steady_state_sweep
export plot_1D_solutions_branch
export follow_branch

# Error hint for extensions stubs
function __init__()
    Base.Experimental.register_error_hint(
        _error_hinter("OrdinaryDiffEq", :TimeEvolution, follow_branch), MethodError
    )
    Base.Experimental.register_error_hint(
        _error_hinter("OrdinaryDiffEq", :TimeEvolution, plot_1D_solutions_branch),
        MethodError,
    )
    Base.Experimental.register_error_hint(
        _error_hinter("SteadyStateDiffEq", :SteadyStateDiffEqExt, steady_state_sweep),
        MethodError,
    )
    for func in [
        plot_spaghetti,
        plot_eigenvalues,
        plot_rotframe_jacobian_response,
        plot_phase_diagram,
        plot_linear_response,
    ]
        Base.Experimental.register_error_hint(
            _error_hinter("Plots", :PlotsExt, func), MethodError
        )
    end
    Base.Experimental.register_error_hint(
        _error_hinter("HarmonicBalance", :HarmonicBalanceExt, get_linear_response),
        MethodError,
    )
    return nothing
end

end
