using HarmonicSteadyState, QuantumCumulants

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
        χ = get_forward_transmission_response(result, Ω_range, 3)

        κ_ext = 0.05
        S21 = 1 .- χ * κ_ext / 2

        @test minimum(abs.(S21)) > 0
        @test maximum(abs.(S21)) < 1
    end
end
