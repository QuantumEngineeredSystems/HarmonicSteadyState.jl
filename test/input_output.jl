using HarmonicSteadyState, QuantumCumulants, Test

@testset "magnon-polariton" begin
    # Hilbertspace
    hm = FockSpace(:magnon)
    hc = FockSpace(:polariton)
    h = hm ⊗ hc

    # Operators
    @qnumbers m::Destroy(h, 1) c::Destroy(h, 2)
    @rnumbers Δ Vk Ωd γm γk
    param = [Δ, Vk, Ωd, γm, γk]

    H_RWA_sym = (
        Δ * m' * m +
        Δ / 2 * c' * c +
        Vk * m * c' * c' +
        Vk * m' * c * c +
        (Ωd * m + Ωd * m')
    )
    ops = [m, m', c, c']
    eqs_RWA = complete(meanfield(ops, H_RWA_sym, [m, c]; rates=[γm, γk], order=1))

    harmonic_eq = HarmonicEquation(eqs_RWA, param)

    fixed = (Δ => 0, Vk => 0.0002, γm => 0.1, γk => 0.01)
    drive_range = range(0, 1, 100)
    varied = (Ωd => drive_range)
    result = get_steady_states(harmonic_eq, TotalDegree(), varied, fixed)

    Ω_range = range(-0.2, 0.2, 500)

    @testset "S21" begin
        m = result.problem.variables[1]
        χ = get_susceptibility(result, 1, Ω_range, 3)

        κ_ext = 0.05
        S21 = 1 .- χ * κ_ext / 2

        S21_test = get_forward_transmission_response(
            result, 1, Ω_range, 3, κ_ext; class="stable"
        )

        @test minimum(abs.(S21)) > 0
        @test maximum(abs.(S21)) < 1

        @test all(iszero, real(S21[1:250, :] - reverse(S21[251:end, :]; dims=1)))
        @test all(iszero, imag(S21[1:250, :] + reverse(S21[251:end, :]; dims=1)))

        @test all(iszero, S21_test - S21)

        @testset "peaks" begin
            using Peaks
            absline = -1 .* abs.(S21[:,end])
            idxs, _ = findmaxima(absline)
            @test length(idxs) == 2
        end
    end

    @testset "make_S transformation" begin
        using HarmonicSteadyState.LinearResponse: make_S
        using LinearAlgebra: I
        # Test S matrix properties
        N = 4  # Even number
        S, Sinv = make_S(N)

        @test size(S) == (N, N)
        @test size(Sinv) == (N, N)
        @test S * Sinv ≈ I(N) rtol = 1e-12  # Should be inverse
        @test Sinv * S ≈ I(N) rtol = 1e-12

        # Test assertion for odd N
        @test_throws AssertionError make_S(3)
    end
end
