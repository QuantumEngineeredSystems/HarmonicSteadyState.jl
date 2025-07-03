@testset "Concretely typed" begin
    CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
    julia_version = VERSION >= v"1.11.0-DEV.0" # fails on 1.11
    if !CI && !julia_version
        using HarmonicSteadyState

        using CheckConcreteStructs

        all_concrete(HarmonicSteadyState.WarmUp)
        all_concrete(HarmonicSteadyState.TotalDegree)
        all_concrete(HarmonicSteadyState.Polyhedral)
        all_concrete(HarmonicSteadyState.Result)
        all_concrete(HarmonicSteadyState.HomotopyContinuationProblem)
        all_concrete(HarmonicSteadyState.HarmonicEquation)
        all_concrete(HarmonicSteadyState.AdiabaticSweep)

        all_concrete(HarmonicSteadyState.LinearResponse.Lorentzian)
        all_concrete(HarmonicSteadyState.LinearResponse.ResponseMatrix)
        all_concrete(HarmonicSteadyState.LinearResponse.JacobianSpectrum)
    end
end

@testset "ExplicitImports" begin
    using ExplicitImports

    @test check_no_implicit_imports(HarmonicSteadyState) == nothing
    @test check_all_explicit_imports_via_owners(HarmonicSteadyState) == nothing
    # @test check_all_explicit_imports_are_public(HarmonicSteadyState) == nothing
    @test check_no_stale_explicit_imports(HarmonicSteadyState) == nothing
    @test check_all_qualified_accesses_via_owners(HarmonicSteadyState) == nothing
    # @test check_all_qualified_accesses_are_public(HarmonicSteadyState) == nothing
    @test check_no_self_qualified_accesses(HarmonicSteadyState) == nothing
end

@testset "Code quality" begin
    using Aqua
    using OrdinaryDiffEqTsit5, SteadyStateDiffEq, Plots, HarmonicBalance

    TimeEvolution = Base.get_extension(HarmonicSteadyState, :TimeEvolution)
    SteadyStateDiffEqExt = Base.get_extension(HarmonicSteadyState, :SteadyStateDiffEqExt)
    PlotsExt = Base.get_extension(HarmonicSteadyState, :PlotsExt)
    HarmonicBalanceExt = Base.get_extension(HarmonicSteadyState, :HarmonicBalanceExt)

    Aqua.test_ambiguities([HarmonicSteadyState])
    using HarmonicSteadyState.HomotopyContinuation: ModelKit
    Aqua.test_all(
        HarmonicSteadyState;
        piracies=(treat_as_own=[ModelKit.Variable, ModelKit.System],),
        ambiguities=false,
        persistent_tasks=false,
    )
    for mod in [TimeEvolution, SteadyStateDiffEqExt, PlotsExt, HarmonicBalanceExt]
        @test check_no_stale_explicit_imports(mod) == nothing
        @test check_all_explicit_imports_via_owners(mod) == nothing
        # Aqua.test_ambiguities(mod)
        Aqua.test_all(
            mod;
            deps_compat=false,
            ambiguities=false,
            piracies=false,
            stale_deps=false,
            project_extras=false,
            persistent_tasks=false,
        )
    end
end

if VERSION < v"1.12.0-beta"
    @testset "Code linting" begin
        using JET
        JET.test_package(HarmonicSteadyState; target_defined_modules=true)
    end
end
