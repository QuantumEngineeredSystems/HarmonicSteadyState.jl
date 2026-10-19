using HarmonicBalance, HarmonicSteadyState, LinearAlgebra, Test
HSS = HarmonicSteadyState

# x'' + γx' + x³ = F cos(ωt) with an [ω, 3ω] ansatz. The references are the Floquet exponents
# of the periodic orbits of this equation itself, found by shooting on it: at ω = 0.8 its
# only orbit is stable with both exponents at -γ/2 = -0.025, and at ω = 1.36 the middle of
# its three orbits is a saddle with exponent 0.323.
@variables α ω ω0 F γ η t x(t) y(t)

function duffing(harmonics; degree)
    eq = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
    )
    add_harmonic!(eq, x, harmonics)
    return get_harmonic_equations(eq; degree)
end
fixed = (α => 1.0, ω0 => 0.0, F => 1.0, γ => 0.05)

function spectra(res, index)
    spectrum = HSS.stability_spectrum(res)
    return [
        spectrum(HSS.get_variable_solutions(res; branch, index)) for
        branch in findall(res.classes["physical"][index])
    ]
end

@testset "slow flow misjudges a stable orbit" begin
    res = get_steady_states(
        duffing([ω, 3ω]; degree=2), ω => [0.8], fixed; show_progress=false
    )
    @test HSS.stability_spectrum(res) isa HSS.FirstOrderSpectrum
    @test count(res.classes["physical"][1]) == 1
    @test !any(res.classes["stable"][1]) # the orbit is stable, see above
end

@testset "second-order spectrum gives the Floquet exponents" begin
    res = get_steady_states(
        duffing([ω, 3ω]; degree=3), ω => [0.8], fixed; show_progress=false
    )
    @test HSS.stability_spectrum(res) isa HSS.SecondOrderSpectrum
    @test count(res.classes["stable"][1]) == 1
    @test !any(res.classes["Hopf"][1])
    λs = only(spectra(res, 1))
    @test length(λs) == 8 # two per harmonic variable
    @test all(λ -> isapprox(real(λ), -0.025; atol=1e-6), λs)

    res = get_steady_states(
        duffing([ω, 3ω]; degree=3), ω => [1.36], fixed; show_progress=false
    )
    @test count(res.classes["physical"][1]) == 3
    @test count(res.classes["stable"][1]) == 2
    saddle = maximum(λs -> maximum(real, λs), spectra(res, 1))
    @test saddle ≈ 0.323 atol = 0.005
end

@testset "first-order variables leave the pencil singular" begin
    # y relaxes towards x and is first order, so it has no second slow-time derivatives
    eqs = DifferentialEquation(
        [
            d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) + η * y ~ F * cos(ω * t),
            d(y, t) + y ~ x,
        ],
        [x, y],
    )
    add_harmonic!(eqs, x, ω)
    add_harmonic!(eqs, y, ω)
    heq = get_harmonic_equations(eqs; degree=3)
    res = get_steady_states(
        heq,
        ω => [1.2],
        (α => 1.0, ω0 => 1.0, F => 0.1, γ => 0.05, η => 0.01);
        show_progress=false,
    )
    @test HSS.stability_spectrum(res) isa HSS.SecondOrderSpectrum
    # two eigenvalues for each of u1, v1 of x, one for each of u2, v2 of y
    @test all(λs -> length(λs) == 6, spectra(res, 1))
end
