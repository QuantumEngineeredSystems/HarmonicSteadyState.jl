using HarmonicSteadyState
using HarmonicBalance
using Test

using Random
const SEED = 0x8f88209c
Random.seed!(SEED)

# JET runs in its own CI job: JET >= 0.10 requires Julia 1.12, so it cannot be
# installed on the LTS runner. GROUP=Core skips it, GROUP=JET runs only it.
const GROUP = get(ENV, "GROUP", "All")

if GROUP in ("All", "Core")
    @testset "Code quality" begin
        include("code_quality.jl")
    end

    @testset "API" begin
        include("API.jl")
        include("Problem.jl")
    end

    @testset "Computing steady states" begin
        @testset "parametron" begin
            include("steady_states/parametron.jl")
        end
        # include("steady_states/krylov.jl")
        include("steady_states/methods.jl")
    end

    @testset "Processing solutions" begin
        include("Jacobian.jl")
        include("transform_solutions.jl")
        include("sorting.jl")
        include("classification.jl")
    end

    @testset "Plotting" begin
        include("plotting.jl")
    end

    @testset "Linear response" begin
        include("linear_response.jl")
        include("input_output.jl")
    end

    @testset "Limit cycle" begin
        include("limit_cycle.jl")
    end

    @testset "extensions" begin
        @testset "QuantumCumulants extension" begin
            include("extensions/QuantumCumulantsExt.jl")
        end
        @testset "Time evolution extension" begin
            include("extensions/time_evolution.jl")
            include("extensions/hysteresis_sweep.jl")
        end
        @testset "SteadyState extension" begin
            include("extensions/SteadyStateDiffEqExt.jl")
        end
    end
end

if GROUP in ("All", "JET")
    @testset "Code linting" begin
        include("jet.jl")
    end
end
