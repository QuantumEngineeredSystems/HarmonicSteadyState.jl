using QuantumCumulants, HarmonicSteadyState, Symbolics
# using Plots
using Test

@testset "KPO" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables Δ::Real U::Real F::Real G::Real κ::Real

    H_RWA = (-Δ + U) * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2
    ops = [a, a']

    eqs_RWA = meanfield(ops, H_RWA, [a]; rates=[κ], order=1)
    eqs_completed_RWA = complete(eqs_RWA)

    fixed = (U => 0.001, κ => 0.00, Δ => 0.0)
    varied = (G => range(0.01, 0.02, 10))

    @testset "HomotopyContinuationProblem" begin
        problem = HarmonicSteadyState.HomotopyContinuationProblem(
            eqs_completed_RWA, [Δ, U, G, κ], varied, fixed
        )
        result = get_steady_states(problem, TotalDegree())
        @test sum(all.(get_class(result, "stable"))) == 2

        fixed = (U => 0.001, κ => 0.00, G => 0.01)
        varied = (Δ => range(-0.03, 0.03, 10))
        problem = HarmonicSteadyState.HomotopyContinuationProblem(
            eqs_completed_RWA, [Δ, U, G, κ], varied, fixed
        )
        result = get_steady_states(problem, TotalDegree())
        @test sum(any.(get_class(result, "stable"))) == 3
    end

    @testset "HarmonicEquation" begin
        harmonic_eq = HarmonicSteadyState.HarmonicEquation(eqs_completed_RWA, [Δ, U, G, κ])
        result = get_steady_states(harmonic_eq, varied, fixed)
        @test sum(any.(get_class(result, "stable"))) == 3
    end

    @testset "second order cumulant" begin
        eqs_RWA = meanfield(ops, H_RWA, [a]; rates=[κ], order=2)
        eqs_completed_RWA = complete(eqs_RWA)

        fixed = (U => 0.001, κ => 0.002, G => 0.01)
        varied = (Δ => range(-0.03, 0.01, 100))
        problem_c2 = HarmonicSteadyState.HomotopyContinuationProblem(
            eqs_completed_RWA, [Δ, U, G, κ], varied, fixed
        )
        @test length(problem_c2.variables) == 5

        result = get_steady_states(problem_c2, TotalDegree())
        @test sum(any.(get_class(result, "stable"))) == 5
        classify_solutions!(result, "a⁺aᵣ < 0", "neg photon number")
        @test maximum(
            phase_diagram(result; not_class="neg photon number", class="stable")
        ) == 3
    end
end

@testset "work with rnumbers and cumber" begin
    @testset "@cnumbers" begin
        h = FockSpace(:cavity)
        @qnumbers a::Destroy(h)
        @cnumbers Δ U G κ
        param = [Δ, U, G, κ]

        H_RWA = -Δ * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2
        ops = [a, a']

        eqs = meanfield(ops, H_RWA, [a]; rates=[κ], order=1)

        fixed = (U => 0.001, κ => 0.002)
        varied = (Δ => range(-0.03, 0.03, 10), G => range(1e-5, 0.02, 10))
        problem_c1 = HarmonicSteadyState.HomotopyContinuationProblem(
            complete(eqs), param, varied, fixed
        )
    end
    @testset "@rnumbers" begin
        h = FockSpace(:cavity)
        @qnumbers a::Destroy(h)
        @rnumbers Δ U G κ
        param = [Δ, U, G, κ]

        H_RWA = -Δ * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2
        ops = [a, a']

        eqs = meanfield(ops, H_RWA, [a]; rates=[κ], order=1)

        fixed = (U => 0.001, κ => 0.002)
        varied = (Δ => range(-0.03, 0.03, 50), G => range(1e-5, 0.02, 50))
        problem_c1 = HarmonicSteadyState.HomotopyContinuationProblem(
            complete(eqs), param, varied, fixed
        )
    end
end

@testset "Multiple modes" begin
    # Hilbertspace
    hc = FockSpace(:cavity)
    hm = FockSpace(:motion)
    h = hc ⊗ hm

    # Operators
    @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

    # Parameters
    @rnumbers Δ K F κ ωm g0 Γm

    param = [Δ, K, F, κ, ωm, g0, Γm]

    H_RWA =
        -Δ * a' * a +
        ωm * b' * b +
        K / 2 * (a'^2 * a^2) +
        F * (a' + a) +
        g0 * a' * a * (b + b')
    ops = [a, a', b, b']

    eqs_RWA = meanfield(ops, H_RWA, [a, b]; rates=[κ, Γm], order=1)
    eqs_completed_RWA = complete(eqs_RWA)

    fixed = (K => -0.001, κ => 0.002, ωm => 0.01, g0 => 0.001, Γm => 0.0001)
    varied = (Δ => range(-0.03, 0.03, 100), F => range(1e-5, 0.02, 100))
    problem = HarmonicSteadyState.HomotopyContinuationProblem(
        eqs_completed_RWA, param, varied, fixed
    )
end
