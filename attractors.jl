using HarmonicBalance, HarmonicSteadyState
using Symbolics, RuntimeGeneratedFunctions, QuestBase, FunctionWrappers, SciMLBase
using ModelingToolkit
import FunctionWrappers: FunctionWrapper
using Attractors

@variables γ λ α ω0 ω t
@variables x(t)

nat_eq = d(d(x, t), t) + γ * d(x, t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * x + α * x^3
dEOM = DifferentialEquation(nat_eq, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM);

struct AttractorsProblem{
    ParType<:Number,
    Jac<:JacobianFunction(ComplexF64),
} <: SteadyStateProblem
    variables::Vector{Num}
    parameters::Vector{Num}
    swept_parameters::OrderedDict{Num,Vector{ParType}}
    fixed_parameters::OrderedDict{Num,ParType}
    phase_space::NTuple{N,Vector{ParType}}
    coupled_ode::DynamicalSystemBase.CoupledODEs
    jacobian::Jac
    eom::HarmonicEquation
end

ODEProblem(
    harmonic_eq,
    zeros(2),
    (0, Inf),
    Dict(ω0 => 1.0, γ => 0.1, λ => 0.05, α => 0.01),
)

u_range = range(-0.2, 0.2; length=75)
u_range |> typeof
grid = (
    u_range, # u1
    u_range, # v1
    u_range, # u2
    u_range, # v2
    u_range, # u3
    u_range, # v3
)
mapper = AttractorsViaRecurrences(
    ds,
    grid;
    consecutive_recurrences=5000,
    attractor_locate_steps=5000,
    consecutive_lost_steps=500,
)

fs = basins_fractions(mapper, sampler)
attractors = extract_attractors(mapper)

plot_attractors(attractors; access=SVector(1, 2))

continuation_range = range(1e-4, 0.002, 100)

sampler, = statespace_sampler(grid)
algo = AttractorSeedContinueMatch(mapper)

reinit!(ds)
fractions_cont, attractors_cont = global_continuation(
    algo, continuation_range, :F₂, sampler; samples_per_parameter=50
)


# eqs = rearrange_standard(harmonic_eq).equations
# rhs = Num[eq.lhs for eq in eqs]
# vars = get_variables(harmonic_eq)

# function compile_function(
#     rhs::Vector{Num}, variables::Vector{Num}; rules=Dict()
# )::RuntimeGeneratedFunction
#     rhss = QuestBase.substitute_all.(rhs, Ref(rules)) # Ref makes sure only mat is broadcasted
#     rhsf = Symbolics.build_function(rhss, variables; expression=Val(false))
#     return rhsf isa Tuple ? first(rhsf) : rhsf
# end

# rgf = compile_function(
#     rhs, vcat(vars, ω); rules=Dict(ω0 => 1.0, γ => 0.1, λ => 0.05, α => 0.01)
# ) # RuntimeGeneratedFunction

# function f(u,p,t)
#     rgf(vcat(u,p))
# end
# ff = FunctionWrapper{Vector{Float64},Tuple{Vector{Float64},Vector{Float64},Float64}}(f)
# ff(rand(2), rand(1), 0.0)

# CoupledODEs(ff, rand(2) , rand(1))
# SciMLBase.isinplace(ff, 100)
