# Attractors on an averaged parametron

Find basin fractions for the two stable fixed points of a parametrically driven Duffing oscillator's single-harmonic slow flow. HarmonicBalance derives and rearranges the averaged equations; Symbolics compiles the numerical right-hand side once. Attractors maps prescribed, seeded initial conditions with its recurrence algorithm. An independent HarmonicSteadyState solve supplies the stable fixed points used for correctness checks.

This is an experiment, not a new `SteadyStateProblem` subtype or supported solver integration. Attractors is a dependency only of this environment. The physics here is the averaged harmonic approximation, not the exact driven oscillator dynamics.

Run from the repository root on Julia 1.13:

```sh
julia --project=experiments/attractors -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/attractors experiments/attractors/run.jl
```

The experiment uses the checkout through `[sources]`; HarmonicBalance and QuestBase follow their default Git branches. Preserve the ignored generated Manifest and record `Pkg.status(; mode=Pkg.PKGMODE_MANIFEST)` when sharing measurements.

The dimensionless model parameters are `ω=ω₀=1`, `γ=0.1`, `λ=0.3`, and `α=0.1`. Initial conditions are 64 uniform samples in `[-2,2]²` from `Xoshiro(0x8f88209c)`. A `101×101` recurrence grid, transient time 200, step interval 1, 100 consecutive recurrences, and 100 attractor-location steps are explicit numerical choices. The integrator is Tsit5 with absolute and relative tolerances `1e-9`; unresolved labels cause the correctness checks to fail.

The runner matches attractors by coordinates, not label numbers, and compares them with stable solutions. It checks the reference residuals, complete matching at these parameters, and normalized sample fractions. Grid-based matching tolerance is three grid spacings. Fractions estimate the sampled region and do not establish global basin volumes or absence of unsampled attractors. Repeat with `AttractorsExperiment.run(; count=128, resolution=201)` to assess sampling and grid sensitivity.

Benchmarks reconstruct all mutable integration and recurrence state for each warmed sample. Their times exclude symbolic derivation, compilation, and the reference steady-state solve. CI publishes results and the resolved Manifest as artifacts; timing has no pass/fail threshold.

The runner uses the [Attractors 2 basin-map and prescribed-initial-condition API](https://juliadynamics.github.io/Attractors.jl/dev/api/), without plotting dependencies or global random-state mutation for sampling.
