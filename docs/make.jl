CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using HarmonicSteadyState
using Documenter
using Plots, SteadyStateDiffEq, OrdinaryDiffEqTsit5

TimeEvolution = Base.get_extension(HarmonicSteadyState, :TimeEvolution)
SteadyStateDiffEqExt = Base.get_extension(HarmonicSteadyState, :SteadyStateDiffEqExt)
PlotsExt = Base.get_extension(HarmonicSteadyState, :PlotsExt)

include("pages.jl")

makedocs(;
    sitename="HarmonicSteadyState.jl",
    authors="Quest group",
    modules=[
        HarmonicSteadyState,
        TimeEvolution,
        SteadyStateDiffEqExt,
        HarmonicSteadyState.LinearResponse,
        PlotsExt,
    ],
    format=Documenter.HTML(;
        canonical="https://quantumengineeredsystems.github.io/HarmonicSteadyState.jl/stable/",
    ),
    pages=pages,
    clean=true,
    linkcheck=true,
    warnonly=:missing_docs,
    draft=(!CI),
    doctest=false,  # We test it in the CI, no need to run it here
)

if CI
    deploydocs(;
        repo="github.com/QuantumEngineeredSystems/HarmonicSteadyState.jl",
        devbranch="main",
        target="build",
        branch="gh-pages",
        push_preview=true,
    )
end
