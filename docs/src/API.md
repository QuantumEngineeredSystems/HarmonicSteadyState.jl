```@meta
CollapsedDocStrings = true
```

# API

```@contents
Pages = ["API.md"]
Depth = 2:3
```

## Finding steady states

```@docs
get_steady_states
HarmonicSteadyState.SteadyStateProblem
HarmonicSteadyState.SteadyStateMethods
HarmonicSteadyState.Result
```

### Homotopy Continuation

```@docs
HarmonicSteadyState.HomotopyContinuationProblem
HarmonicSteadyState.HomotopyContinuationMethods
WarmUp
TotalDegree
Polyhedral
```

## Analyze solutions

### Access steady states

```@docs
get_solutions
get_single_solution
attractors
transform_solutions
```

### Classify steady states

```@docs
classify_solutions!
get_class
filter_result!
phase_diagram
```

### Steady state plotting

To use these plotting functions, you need to have the `Plots` package installed in the same enverionment and loaded.

```@docs
plot(::HarmonicSteadyState.Result, args...; kwargs...)
plot!(::HarmonicSteadyState.Result, args...; kwargs...)
plot_phase_diagram
plot_spaghetti
```

## Linear response

```@autodocs
Modules = [HarmonicSteadyState.LinearResponse]
Private = false
Order = [:function]
```

### Linear response plotting

```@docs
plot_eigenvalues
plot_linear_response
plot_rotframe_jacobian_response
```

## Limit-cycle methods

```@docs
get_limit_cycles
add_pairs!
get_cycle_variables
```

### OrdinaryDiffEq

```@docs
AdiabaticSweep
follow_branch
plot_1D_solutions_branch
```

```@autodocs; canonical=false
Modules = [Base.get_extension(HarmonicSteadyState, :TimeEvolution)]
Private = false
Order = [:function]
```

### SteadyStateSweep

```@docs
steady_state_sweep
```
