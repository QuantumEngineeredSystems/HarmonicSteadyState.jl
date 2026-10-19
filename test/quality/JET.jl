using HarmonicSteadyState
using Test
using JET

@testset "JET" begin
    JET.test_package(HarmonicSteadyState; target_modules=(HarmonicSteadyState,))
end
