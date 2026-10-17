using HarmonicSteadyState
using OrdinaryDiffEqTsit5
using Test
using Symbolics: Symbolics, Num

rhs_allocations(f, du, u, time) = @allocated f(du, u, nothing, time)

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
    @test du != expected[2]

    @variables u1 v1
    flat_initial = OrderedDict(fixed...)
    flat_initial[ω] = 1.0
    flat_initial[v1] = -0.3
    flat_initial[u1] = 0.2
    @test ODEProblem(harmonic_eq, flat_initial; timespan=tspan).u0 == state

    # Multiple independent sweeps must follow their matching parameter order.
    combined = sweep + AdiabaticSweep(F => (1e-3, 2e-3), tspan)
    multiple = ODEProblem(harmonic_eq, fixed; sweep=combined, u0=state, timespan=tspan)
    standard = rearrange_standard(harmonic_eq)
    for time in (0.0, 5.0, 10.0)
        rules = Dict{Num,Any}(fixed)
        merge!(rules, Dict(zip(variables, state)))
        for parameter in keys(combined)
            rules[parameter] = combined[parameter](time)
        end
        reference = [
            real(Symbolics.value(Symbolics.substitute(eq.lhs, rules; fold=Val(true)))) for
            eq in standard.equations
        ]
        multiple.f(du, state, nothing, time)
        @test du ≈ reference
    end
end

@testset "finite-time stability around zero" begin
    using OrderedCollections: OrderedDict
    using Symbolics: Symbolics
    @variables τ z(τ) rate frequency
    equation = DifferentialEquation(d(z, τ, 2) + rate * d(z, τ) + frequency^2 * z ~ 0, z)
    add_harmonic!(equation, z, frequency)
    harmonics = get_harmonic_equations(equation)
    variables = Symbolics.get_variables(harmonics)
    values = OrderedDict(rate => 0.2, frequency => 1.0)
    for variable in variables
        values[variable] = 0.0
    end
    @test HarmonicSteadyState.is_stable(values, harmonics; timespan=(0.0, 100.0))
    values[rate] = -0.2
    @test !HarmonicSteadyState.is_stable(values, harmonics; timespan=(0.0, 100.0))
end

@testset "stability checks independent perturbation directions" begin
    using OrderedCollections: OrderedDict
    using Symbolics: Symbolics
    @variables τ q(τ) r(τ)
    equation = DifferentialEquation([d(q, τ) ~ -r, d(r, τ) ~ -q], [q, r])
    add_harmonic!(equation, q, 0)
    add_harmonic!(equation, r, 0)
    harmonics = get_harmonic_equations(equation)
    values = OrderedDict(variable => 0.0 for variable in Symbolics.get_variables(harmonics))
    # The uniform direction is stable, but the orthogonal direction is unstable.
    @test !HarmonicSteadyState.is_stable(values, harmonics; timespan=(0.0, 10.0))
end

@testset "stability checks are independent of the amplitude scale" begin
    using OrderedCollections: OrderedDict
    using Symbolics: Symbolics
    @variables τ q(τ) width
    # Bistable flow: zero is unstable and relaxes onto the stable states ±width.
    # With width far below one, both lie within any fixed absolute tolerance of zero.
    equation = DifferentialEquation(d(q, τ) ~ q - q^3 / width^2, q)
    add_harmonic!(equation, q, 0)
    harmonics = get_harmonic_equations(equation)
    variable = only(Symbolics.get_variables(harmonics))
    for (state, stable) in ((0.0, false), (1e-2, true), (-1e-2, true))
        values = OrderedDict(width => 1e-2, variable => state)
        @test HarmonicSteadyState.is_stable(values, harmonics; timespan=(0.0, 20.0)) ==
            stable
    end
end
