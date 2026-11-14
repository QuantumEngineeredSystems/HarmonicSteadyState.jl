module AttractorsExperiment

using Attractors:
    Attractors,
    CoupledODEs,
    BasinMapRecurrences,
    PrescribedICs,
    basins_fractions_labels,
    extract_attractors
using BenchmarkTools: @benchmark, median
using HarmonicBalance: HarmonicBalance
using HarmonicSteadyState: HarmonicSteadyState, TotalDegree, get_steady_states
using LinearAlgebra: norm
using OrdinaryDiffEqTsit5: Tsit5
using Random: Xoshiro, rand
using Symbolics: Symbolics, @variables
using Test: @test, @testset

const SEED = UInt32(0x8f88209c)

function model()
    @variables γ λ α ω0 ω t x(t)
    natural = HarmonicBalance.DifferentialEquation(
        HarmonicBalance.d(x, t, 2) +
        γ * HarmonicBalance.d(x, t) +
        ω0^2 * (1 - λ * cos(2 * ω * t)) * x +
        α * x^3 ~ 0,
        x,
    )
    HarmonicBalance.add_harmonic!(natural, x, ω)
    harmonics = HarmonicBalance.get_harmonic_equations(natural)
    fixed = Dict(γ => 0.1, λ => 0.3, α => 0.1, ω0 => 1.0)
    steady = get_steady_states(
        harmonics,
        TotalDegree(; seed=SEED, thread=false),
        Dict(ω => [1.0]),
        fixed;
        show_progress=false,
    )
    stable = [real.(root) for root in values(only(HarmonicSteadyState.attractors(steady)))]

    slow_flow = HarmonicBalance.rearrange_standard(harmonics)
    variables = Symbolics.get_variables(slow_flow)
    rules = merge(fixed, Dict(ω => 1.0))
    rhs = [Symbolics.substitute(eq.lhs, rules) for eq in slow_flow.equations]
    _, compiled! = Symbolics.build_function(rhs, variables; expression=Val(false))
    return compiled!, stable
end

function map_basins(compiled!, initial_conditions, grid)
    # Construct fresh mutable integration and recurrence state on every run.
    f!(du, u, p, t) = compiled!(du, u)
    system = CoupledODEs(
        f!,
        first(initial_conditions),
        nothing;
        diffeq=(alg=Tsit5(), abstol=1e-9, reltol=1e-9),
    )
    mapper = BasinMapRecurrences(
        system,
        grid;
        Δt=1.0,
        Ttr=200.0,
        consecutive_recurrences=100,
        attractor_locate_steps=100,
        consecutive_lost_steps=100,
        maximum_iterations=10_000,
    )
    fractions, labels = basins_fractions_labels(
        mapper, PrescribedICs(initial_conditions); show_progress=false
    )
    return (; fractions, labels, attractors=extract_attractors(mapper))
end

function run(; count=64, resolution=101, samples=20)
    compiled!, stable = model()
    rng = Xoshiro(SEED)
    initial_conditions = [4 .* rand(rng, 2) .- 2 for _ in 1:count]
    grid = ntuple(_ -> range(-2.0, 2.0; length=resolution), 2)
    result = map_basins(compiled!, initial_conditions, grid)
    # Match point attractors by coordinates, never by package-specific labels.
    centers = [sum(points) / length(points) for points in values(result.attractors)]
    tolerance = 3 * step(first(grid))
    @testset "recurrence attractors agree with stable steady states" begin
        @test length(stable) == 2
        @test length(centers) == length(stable)
        @test all(>(0), result.labels)
        @test sum(values(result.fractions)) ≈ 1
        @test all(
            center -> minimum(norm(center-root) for root in stable) < tolerance, centers
        )
        @test all(
            root -> minimum(norm(root-center) for center in centers) < tolerance, stable
        )
        @test all(values(result.attractors)) do points
            all(point -> minimum(norm(point-root) for root in stable) < tolerance, points)
        end
        for root in stable
            residual = similar(root)
            compiled!(residual, root)
            @test norm(residual, Inf) < 1e-8
        end
    end
    println(
        "Julia ",
        VERSION,
        "; threads=",
        Threads.nthreads(),
        "; Attractors=",
        pkgversion(Attractors),
        "; seed=",
        SEED,
    )
    println("Samples=", count, "; region=[-2,2]²; grid=", resolution, "²; Ttr=200; Δt=1")
    println("Stable steady states: ", stable)
    println("Recurrence centers: ", centers)
    println("Basin fractions (sample estimates): ", result.fractions)
    trial = @benchmark map_basins($compiled!, $initial_conditions, $grid) samples=samples evals=1 seconds=60
    estimate = median(trial)
    println(
        "Fresh basin-map run: median ms=",
        round(estimate.time/1e6; digits=3),
        "; bytes=",
        estimate.memory,
        "; allocations=",
        estimate.allocs,
    )
    println(
        "Timing excludes harmonic derivation, compilation, and the steady-state reference solve.",
    )
    return result
end

end

AttractorsExperiment.run()
