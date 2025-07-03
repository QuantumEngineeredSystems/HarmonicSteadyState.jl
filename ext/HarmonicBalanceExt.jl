module HarmonicBalanceExt

using DocStringExtensions: TYPEDSIGNATURES
using Symbolics: Symbolics, get_variables, Num
using HarmonicBalance: HarmonicBalance, get_Jacobian
using QuestBase:
    QuestBase,
    DifferentialEquation,
    declare_variable,
    get_independent_variables,
    substitute_all,
    d,
    var_name
using HarmonicSteadyState: HarmonicSteadyState
using HarmonicSteadyState: LinearResponse
using HarmonicSteadyState.LinearResponse: ResponseMatrix
using ProgressMeter: ProgressMeter, Progress, next!

"""
$(TYPEDSIGNATURES)

Obtain the symbolic linear response matrix of a `diff_eq` corresponding to a perturbation frequency `freq`.
This routine cannot accept a `HarmonicEquation` since there, some time-derivatives are already dropped.
`order` denotes the highest differential order to be considered.

"""
function get_response_matrix(diff_eq::DifferentialEquation, freq::Num; order=2)::Matrix
    Symbolics.@variables T
    time = get_independent_variables(diff_eq)[1]

    eom = HarmonicBalance.harmonic_ansatz(diff_eq, time)

    # replace the time-dependence of harmonic variables by slow time BUT do not drop any derivatives
    eom = HarmonicBalance.slow_flow(eom; fast_time=time, slow_time=T, degree=order + 1)

    eom = HarmonicBalance.fourier_transform(eom, time)

    # get the response matrix by summing the orders
    # M = Symbolics.jacobian(eom.equations, get_variables(eom))
    M = get_Jacobian(eom.equations, get_variables(eom))
    for n in 1:order
        # M += (im * freq)^n * Symbolics.jacobian(eom.equations, d(get_variables(eom), T, n))
        M += (im * freq)^n * get_Jacobian(eom.equations, d(get_variables(eom), T, n))
    end
    M = substitute_all(
        M, [var => declare_variable(var_name(var)) for var in get_variables(eom)]
    )
    return M
end

"Get the response matrix corresponding to `res`.
Any substitution rules not specified in `res` can be supplied in `rules`."
function HarmonicSteadyState.LinearResponse.ResponseMatrix(
    res::HarmonicSteadyState.Result; rules=Dict()
)
    eom = source(res.problem)
    if isnothing(eom)
        error("Cannot get the response matrix of the second order natural differential equations of a result with a problem with no source.")
    elseif !isa(source(eom),QuestBase.DifferentialEquation)
        error("Cannot get the response matrix the second order natural differential equations of a result where the source system is not second order differential equation.")
    end

    # get the symbolic response matrix
    Symbolics.@variables Δ
    M = get_response_matrix(source(eom), Num(Δ))
    M = QuestBase.substitute_all(M, merge(res.fixed_parameters, rules))
    symbols = HarmonicSteadyState._free_symbols(res)

    compiled_M = map(M) do el
        args = cat(symbols, [Δ]; dims=1)
        f_re = Symbolics.build_function(el.re, args; expression=Val{false})
        f_im = Symbolics.build_function(el.im, args; expression=Val{false})
        (args...) -> f_re(args...) + im * f_im(args...)
    end
    return ResponseMatrix(compiled_M, symbols, eom.variables)
end

"""
$(TYPEDSIGNATURES)

Calculate the linear response of the system for a given branch.
Evaluates the linear response by solving the linear response ODE for each stable solution
and input frequency in the given range.

# Arguments
- `res`: Result object containing the system's solutions
- `nat_var::Num`: Natural variable to evaluate in the response
- `Ω_range`: Range of frequencies to evaluate
- `branch::Int`: Branch number to analyze
- `show_progress=true`: Whether to show a progress bar

# Returns
- Array{P,2}: Response matrix where rows correspond to frequencies and columns to stable solutions
"""
function HarmonicSteadyState.get_linear_response(
    res::HarmonicSteadyState.Result{D,S,P},
    nat_var::Num,
    Ω_range,
    branch::Int;
    show_progress=true,
) where {D,S,P}
    stable = HarmonicSteadyState.get_class(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    response = ResponseMatrix(res) # the symbolic response matrix
    C = Array{P,2}(undef, length(Ω_range), sum(stable))

    # note: this could be optimized by not grabbing the entire huge dictionary every time
    if show_progress
        bar = Progress(
            length(C);
            dt=1,
            desc="Solving the linear response ODE for each solution and input frequency ... ",
            barlen=50,
        )
    end
    for j in findall(stable)

        # get response for each individual point
        s = HarmonicSteadyState.get_single_solution(res; branch=branch, index=j)
        for i in 1:(size(C)[1])
            C[i, j] = HarmonicSteadyState.LinearResponse.get_response(
                response, s, Ω_range[i]
            )
        end
        show_progress ? next!(bar) : nothing
    end
    return C
end

end # module
