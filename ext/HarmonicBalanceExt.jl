module HarmonicBalanceExt

using DocStringExtensions
using Symbolics: Symbolics, get_variables, Num
using HarmonicBalance: HarmonicBalance, get_Jacobian
using QuestBase:
    DifferentialEquation,
    declare_variable,
    get_independent_variables,
    substitute_all,
    d,
    var_name
using HarmonicSteadyState: HarmonicSteadyState

"""
$(TYPEDSIGNATURES)

Obtain the symbolic linear response matrix of a `diff_eq` corresponding to a perturbation frequency `freq`.
This routine cannot accept a `HarmonicEquation` since there, some time-derivatives are already dropped.
`order` denotes the highest differential order to be considered.

"""
function HarmonicSteadyState.get_response_matrix(
    diff_eq::DifferentialEquation, freq::Num; order=2
)::Matrix
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

end
