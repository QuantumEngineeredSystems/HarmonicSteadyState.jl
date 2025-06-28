using HarmonicBalance, ModelingToolkit
using Attractors


@variables γ λ α ω0 ω t
@variables x(t)

nat_eq = d(d(x, t), t) + γ * d(x, t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * x + α * x^3
dEOM = DifferentialEquation(nat_eq, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM);

using OrdinaryDiffEq
param_1d = Dict(
    Δ => Delta / omega0, γ => gamma / omega0, F₁ => F1_fix, δ => δ_exp / omega0, F₂ => 1e-4)
prob = ODEProblem(harmonic_eq, SA[rand(6)...] ./ 4, (0.0, 100.0), param_1d; in_place=false)
diffeq = (alg=Vern7(),abstol = 1e-6,reltol = 1e-6)
ds = CoupledODEs(prob, #=diffeq=#)

set_parameter!(ds, :F₂, 0.0007)

u_range = range(-0.2, 0.2; length=75)
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
