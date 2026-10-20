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
    @test length(λs) == 2 # one per Floquet exponent of x
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
    # two Floquet exponents of x, one of y
    @test all(λs -> length(λs) == 3, spectra(res, 1))
end

@testset "copies of an exponent past the highest harmonic are discarded" begin
    # x'' + 0.005x' + x + x³ = 0.0025 cos(ωt) with [ω, 2ω, 3ω]: near the fold of the upper
    # branch, a copy of the stable exponent -γ/2 shifted by 2ω turns unstable in the
    # truncated quadratic problem; shooting confirms the orbit is stable.
    eq = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
    )
    add_harmonic!(eq, x, [ω, 2ω, 3ω])
    res = get_steady_states(
        get_harmonic_equations(eq; degree=3),
        TotalDegree(),
        ω => [1.0779],
        (α => 1.0, ω0 => 1.0, F => 0.0025, γ => 0.005);
        show_progress=false,
    )
    spectrum = HSS.stability_spectrum(res)
    upper = argmax(
        b -> abs(real(res.solutions[1][b][2])), findall(res.classes["physical"][1])
    )
    s = HSS.get_variable_solutions(res; branch=upper, index=1)
    λs = filter(isfinite, eigvals(HSS.companion_pencil(spectrum, s)...))
    @test maximum(real, λs) > 0 # the spurious copy
    @test maximum(real, spectrum(s)) ≈ -0.0025 atol = 1e-8
    @test res.classes["stable"][1][upper]

    # [ω, 3ω] at ω = 0.68: the orbit breaks its symmetry through a real exponent 0.139, which
    # an odd ansatz sees as a conjugate pair at the edge Im λ ≈ ±ω of the strip
    res = get_steady_states(
        duffing([ω, 3ω]; degree=3), ω => [0.68], fixed; show_progress=false
    )
    @test count(res.classes["physical"][1]) == 1
    λs = only(spectra(res, 1))
    @test maximum(real, λs) ≈ 0.114 atol = 0.005
    @test !any(res.classes["stable"][1])
    @test !any(res.classes["Hopf"][1]) # a real exponent, not a Hopf pair
end

@testset "Hopf instability from the second-order spectrum" begin
    # a driven van der Pol oscillator off its locking range sits on a Hopf-unstable focus
    @variables μ
    eq = DifferentialEquation(
        d(x, t, 2) - μ * (1 - x^2) * d(x, t) + ω0^2 * x ~ F * cos(ω * t), x
    )
    add_harmonic!(eq, x, ω)
    res = get_steady_states(
        get_harmonic_equations(eq; degree=3),
        ω => [1.1],
        (μ => 0.1, ω0 => 1.0, F => 0.01);
        show_progress=false,
    )
    @test HSS.stability_spectrum(res) isa HSS.SecondOrderSpectrum
    @test count(res.classes["Hopf"][1]) == count(res.classes["physical"][1]) == 1

    # a saddle has a real unstable exponent, however many complex copies it has
    res = get_steady_states(
        duffing([ω, 3ω]; degree=3), ω => [1.36], fixed; show_progress=false
    )
    @test !any(res.classes["Hopf"][1])
end

@testset "linear response from the second-order spectrum" begin
    using HarmonicSteadyState.LinearResponse: make_S, susceptibility!
    single = (α => 1.0, ω0 => 1.0, F => 0.01, γ => 0.01)
    res = get_steady_states(duffing(ω; degree=3), ω => [0.9], single; show_progress=false)
    @test length(only(eigenvalues(res, 1))) == 2
    @test size(only(eigenvectors(res, 1))) == (2, 2)

    # the Lorentzians peak where the exact response does, which the slow flow misses
    Ω_range = range(0.85, 1.15, 301)
    exact = get_linear_response(res, x, Ω_range, 1; show_progress=false)
    lorentzian = get_jacobian_response(res, x, Ω_range, 1; show_progress=false)
    @test argmax(vec(lorentzian)) == argmax(vec(exact))
    slow = get_steady_states(duffing(ω; degree=2), ω => [0.9], single; show_progress=false)
    slow_flow = get_jacobian_response(slow, x, Ω_range, 1; show_progress=false)
    @test argmax(vec(slow_flow)) != argmax(vec(exact))

    # with J₂ = 0 the susceptibility reduces to the first-order one
    eom = duffing(ω; degree=2)
    slow_spectrum = HSS.stability_spectrum(slow)
    J1 = HSS.compile_matrix(
        HSS._get_J_matrix(eom; order=1),
        HSS._free_symbols(slow);
        rules=slow.fixed_parameters,
    )
    reduced = HSS.SecondOrderSpectrum(
        slow.jacobian, J1, v -> zeros(2, 2), v -> ones(1, 1), 2
    )
    sol = HSS.get_variable_solutions(slow; branch=1, index=1)
    S, Sinv = make_S(2)
    Ωs = range(-0.2, 0.2, 41)
    χ1 = susceptibility!(zeros(ComplexF64, 41), slow_spectrum, sol, 1, Ωs, S, Sinv)
    χ2 = susceptibility!(zeros(ComplexF64, 41), reduced, sol, 1, Ωs, S, Sinv)
    @test χ1 ≈ χ2
end
