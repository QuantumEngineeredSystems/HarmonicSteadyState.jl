using HarmonicSteadyState
using HarmonicBalance
using OrdinaryDiffEqTsit5
using Test
using Symbolics: Symbolics, Num

rhs_allocations(f, du, u, time) = @allocated f(du, u, nothing, time)

"The right-hand side of the rearranged `eom` with `rules` substituted symbolically."
function symbolic_rhs(eom, rules)
    return [
        real(Symbolics.value(Symbolics.substitute(eq.lhs, rules; fold=Val(true)))) for
        eq in rearrange_standard(eom).equations
    ]
end

@variables Ω, γ, λ, F, x, θ, η, α, ω0, ω, t, T, ψ
@variables x(t)

natural_equation =
    d(d(x, t), t) +
    γ * d(x, t) +
    Ω^2 * (1 - λ * cos(2 * ω * t + ψ)) * x +
    α * x^3 +
    η * d(x, t) * x^2
forces = F * cos(ω * t + θ)
dEOM = DifferentialEquation(natural_equation + forces, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM; slow_time=T, fast_time=t);

fixed = (Ω => 1.0, γ => 1e-2, λ => 5e-2, F => 1e-3, α => 1.0, η => 0.3, θ => 0, ψ => 0)

tspan = (0.0, 10)
sweep = AdiabaticSweep(ω => (0.9, 1.1), tspan) # linearly interpolate between two values at two times
ode_problem = ODEProblem(harmonic_eq, fixed; sweep=sweep, u0=[0.01; 0.0], timespan=tspan)
time_soln = solve(ode_problem, Tsit5(); saveat=1);

transform_solutions(time_soln, "sqrt(u1^2+v1^2)", harmonic_eq)

@testset "compiled adiabatic dynamics" begin
    using Symbolics: Symbolics
    using OrderedCollections: OrderedDict

    state = [0.2, -0.3]
    expected = (
        [-0.05426936896251462, -0.026113938441772903],
        [-0.02033407650701693, -0.004546156469083027],
        [0.01014128625550847, 0.014950416193026346],
    )
    for (time, derivative) in zip((0.0, 5.0, 10.0), expected)
        du = zeros(2)
        @test (@inferred ode_problem.f(du, state, nothing, time)) === nothing
        @test du ≈ derivative
    end
    du = zeros(2)
    ode_problem.f(du, state, nothing, 5.0) # warm allocation measurement
    rhs_allocations(ode_problem.f, du, state, 5.0)
    @test rhs_allocations(ode_problem.f, du, state, 5.0) == 0

    # Insertion order must not become state order. Functional state symbols are
    # initial values rather than constant substitutions in the dynamics.
    variables = Symbolics.get_variables(harmonic_eq)
    initial = OrderedDict(fixed...)
    initial[ω] = 1.0
    initial[variables[2]] = -0.3
    initial[variables[1]] = 0.2
    constant_problem = ODEProblem(harmonic_eq, initial; timespan=tspan)
    @test constant_problem.u0 == state
    constant_problem.f(du, [0.4, -0.6], nothing, 5.0)
    @test du ≈ symbolic_rhs(
        harmonic_eq, merge(Dict{Num,Any}(initial), Dict(zip(variables, [0.4, -0.6])))
    )

    @variables u1 v1
    flat_initial = OrderedDict(fixed...)
    flat_initial[ω] = 1.0
    flat_initial[v1] = -0.3
    flat_initial[u1] = 0.2
    @test ODEProblem(harmonic_eq, flat_initial; timespan=tspan).u0 == state

    # Multiple independent sweeps must follow their matching parameter order.
    combined = sweep + AdiabaticSweep(F => (1e-3, 2e-3), tspan)
    multiple = ODEProblem(harmonic_eq, fixed; sweep=combined, u0=state, timespan=tspan)
    for time in (0.0, 5.0, 10.0)
        rules = Dict{Num,Any}(fixed)
        merge!(rules, Dict(zip(variables, state)))
        for parameter in keys(combined)
            rules[parameter] = combined[parameter](time)
        end
        reference = symbolic_rhs(harmonic_eq, rules)
        multiple.f(du, state, nothing, time)
        @test du ≈ reference
    end
end
