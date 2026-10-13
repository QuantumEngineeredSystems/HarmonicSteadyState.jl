# Working with results

A [`Result`](@ref HarmonicSteadyState.Result) stores every computed solution at every point of a parameter
sweep. You can inspect and extract these data directly, without plotting them.
The accessors below also apply to results obtained from harmonic equations or
QuantumCumulants mean-field equations.

## A small example

Consider the steady states of ``\dot u=u-u^3-p`` and ``\dot v=u-v``. Sweep the
parameter ``p`` over three values:

```@example results
using HarmonicSteadyState
using Symbolics: @variables, Num, get_variables

@variables u v p
problem = HarmonicSteadyState.HomotopyContinuationProblem(
    [u - u^3 - p, u - v], [u, v], [p],
    Dict(p => [-0.2, 0.0, 0.2]), Dict{Num,Float64}(),
)
result = get_steady_states(problem, TotalDegree(); show_progress=false)
```

## Inspect the stored data

Julia's `propertynames` lists the available fields. `get_variables` gives the
order of entries in each solution vector, and `size` gives the parameter grid:

```@example results
propertynames(result), get_variables(result), size(result)
```

| Field | Contents |
|:------|:---------|
| `solutions` | An array over the parameter grid. Each entry contains one vector per solution branch. |
| `swept_parameters` | Parameter names and their sweep axes, in grid-dimension order. |
| `fixed_parameters` | Parameter values held fixed throughout the sweep. |
| `classes` | Classification masks such as `"physical"` and `"stable"`, with one flag per branch at each grid point. |
| `problem` | The problem used to compute the result. |
| `jacobian` | The compiled stability Jacobian. |
| `binary_labels` | Integer labels encoding the stable-branch pattern at each grid point. |
| `seed` | The solver seed used for this result. |

For a one-dimensional sweep, `result.solutions[i][branch][variable]` selects a
coordinate. For two dimensions, use
`result.solutions[i, j][branch][variable]`; the first and second indices follow
the first and second keys of `result.swept_parameters`.

Branch identifiers follow the solver's branch sorting. A branch can be absent
at some parameter points; those entries are padded with `NaN`. A computed
solution can also be complex or unstable. Use classification masks when
extracting physical data.

## Extract a point with its parameters

[`get_single_solution`](@ref) returns a dictionary containing all variable,
swept-parameter, and fixed-parameter values for one branch at one grid point:

```@example results
point = get_single_solution(result; branch=1, index=2)
point[u], point[v], point[p]
```

In a two-dimensional sweep, pass `index=(i, j)`. `result[i]` returns the
corresponding dictionaries for every branch at that grid point.

## Extract and filter a branch

[`get_class`](@ref) supplies the stable-point mask for a branch, while
[`get_branches`](@ref) evaluates a named variable or expression along it.
Excluded points remain `NaN`, so the output stays aligned with the sweep axis:

```@example results
axis = swept_parameter(result, p)
stable = get_class(result, 1, "stable")
values = only(get_branches(result, "u"; branches=[1], class="stable"))
collect(zip(axis[stable], values[stable]))
```

Use `class="physical"` to include real solutions whether stable or unstable,
or `class="all"` to inspect every computed solution. Multiple classes select
their intersection; `not_class` excludes additional classes. For complete
solution vectors at each grid point, use [`get_solutions`](@ref):

```@example results
get_solutions(result; class=["physical", "stable"])
```
