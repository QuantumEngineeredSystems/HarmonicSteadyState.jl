"""
$(TYPEDSIGNATURES)

Creates a solution class in `res` using the function `func` (parsed into Symbolics.jl input).
The new class is labeled with `name` and stored under `res.classes[name]`.
By default, only physical (real) solutions are classified, and `false` is returned for the rest.
To also classify complex solutions, set `physical=false`.

# Example
```julia
# solve a previously-defined problem
res = get_steady_states(problem, swept_parameters, fixed_parameters)

# classify, store in result.classes["large_amplitude"]
classify_solutions!(res, "sqrt(u1^2 + v1^2) > 1.0" , "large_amplitude")
```

"""
function classify_solutions!(
    res::Result, func::Union{String,Function}, name::String; physical=true
)
    values = classify_solutions(res, func; physical)
    return res.classes[name] = values
end

"""
$(TYPEDSIGNATURES)

Classifies solutions in `res` using the function `func`.
If `physical` is true, only physical solutions are classified.
"""
function classify_solutions(res::Result, func; physical=true)
    func = isa(func, Function) ? func : _build_substituted(func, res)
    if physical
        f_comp(soln) = _is_physical(soln) && func(real.(soln))
        return transform_solutions(res, f_comp)

    else
        return transform_solutions(res, func)
    end
end

"""
$(TYPEDSIGNATURES)

Returns an array of booleans classifying `branch` in the solutions in `res`
according to `class`.
"""
function get_class(res::Result, branch::Int64, class::String)
    return branch_values = getindex.(res.classes[class], branch)
end

"""
$(TYPEDSIGNATURES)

Returns an array of booleans classifying each branch in the solutions in `res`
according to `class`.
"""
function get_class(soln::Result, class::String)
    return [get_class(soln, b, class) for b in 1:length(first(soln.solutions))]
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` of the Result `res` is physical (real number).
`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_physical(soln::StateDict, res::Result)
    var_values = [getindex(soln, v) for v in res.problem.variables]
    return _is_physical(var_values)
end

_is_physical(soln; im_tol=IM_TOL) = all(x -> !isnan(x) && abs(imag(x)) < im_tol, soln)
_is_physical(res::Result) = classify_solutions(res, _is_physical)

"""
$(TYPEDSIGNATURES)

The spectrum deciding the stability of the steady states in `res`: a
[`SecondOrderSpectrum`](@ref) if the harmonic equations keep the second slow-time
derivatives, a [`FirstOrderSpectrum`](@ref) of `res.jacobian` otherwise. Limit cycles keep
the gauge-fixed Jacobian they were solved with, and problems entered as explicit equations
carry no harmonic equations to take the derivative blocks from.
"""
stability_spectrum(res::Result) = stability_spectrum(res, source(res))
stability_spectrum(res::Result, ::Nothing) = FirstOrderSpectrum(res.jacobian)
function stability_spectrum(res::Result, eom::HarmonicEquation)
    any(hv -> hv.kind == LimitCycleFrequency, eom.variables) &&
        return FirstOrderSpectrum(res.jacobian)
    J2 = _get_J_matrix(eom; order=2)
    all(iszero, J2) && return FirstOrderSpectrum(res.jacobian)
    sym_order = _free_symbols(res)
    rules = res.fixed_parameters
    # each natural variable has as many Floquet exponents as the order of its equation
    natural = [hv.natural_variable for hv in eom.variables]
    exponents = sum(unique(natural)) do x
        return any(!iszero, J2[:, findall(isequal(x), natural)]) ? 2 : 1
    end
    harmonics = Num[hv.kind == ZeroHarmonic ? 0 : hv.ω for hv in eom.variables]
    return SecondOrderSpectrum(
        res.jacobian,
        compile_matrix(_get_J_matrix(eom; order=1), sym_order; rules),
        compile_matrix(J2, sym_order; rules),
        compile_matrix(reshape(harmonics, :, 1), sym_order; rules),
        exponents,
    )
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` of the Result `res` is stable.
Stable solutions are real and have all eigenvalues of their [`stability_spectrum`](@ref)
Re(λ) <= 0.
`im_tol` : an absolute threshold to distinguish real/complex numbers.
`rel_tol`: Re(λ) considered <=0 if real.(λ) < rel_tol*abs(λmax)
"""
function is_stable(soln::StateDict, res::Result; kwargs...)
    return _is_stable(collect(values(soln)), stability_spectrum(res); kwargs...)
end

"""
$(TYPEDSIGNATURES)

Returns a function that checks if a solution is stable.
"""
function _is_stable(res::Result; kwargs...)
    spectrum = stability_spectrum(res)
    return _isit(soln) = _is_stable(soln, spectrum; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` is stable, with `spectrum` its stability spectrum.
"""
function _is_stable(soln, spectrum; rel_tol=1e-10)
    _is_physical(soln) || return false
    λs = spectrum(soln)
    scale = maximum(Iterators.map(abs, λs))
    return all(x -> real(x) < rel_tol * scale, λs)
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` of the Result `res` is Hopf-unstable.
Hopf-unstable solutions are real and unstable, with exactly two eigenvalues of their
[`stability_spectrum`](@ref) in the right half-plane, which are complex conjugates of each other.

`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_Hopf_unstable(soln::StateDict, res::Result)
    return _is_Hopf_unstable(collect(values(soln)), stability_spectrum(res))
end

"""
$(TYPEDSIGNATURES)

Returns a function that checks if a solution is Hopf-unstable.
"""
function _is_Hopf_unstable(res::Result)
    spectrum = stability_spectrum(res)
    return _isit(soln) = _is_Hopf_unstable(soln, spectrum)
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` is Hopf-unstable, with `spectrum` its stability spectrum.
The instability is judged as for the `"stable"` class, so that no solution is both.
"""
function _is_Hopf_unstable(soln, spectrum; rel_tol=1e-10)
    _is_physical(soln) || return false  # the solution is unphysical anyway
    λs = spectrum(soln)
    scale = maximum(Iterators.map(abs, λs))
    unstable = filter(x -> real(x) >= rel_tol * scale, λs)
    return length(unstable) == 2 && abs(conj(unstable[1]) - unstable[2]) < IM_TOL
end

"""
$(TYPEDSIGNATURES)

Create binary classification of the solutions, where each solution point receives an identifier based
on its permutation of stable branches (allows to distinguish between different phases, which may
have the same number of stable solutions). It works by converting each bitstring
`[is_stable(solution_1), is_stable(solution_2), ...,]` into unique labels.
"""
function classify_binaries!(res::Result)
    # mapping of binary string (stable or not) for each solution set to integer.
    # Sensitive to ordering!
    bin_label = bitarr_to_int.(clean_bitstrings(res))

    # renormalize labels with numbers from 1 to length(unique(label))
    for (idx, el) in enumerate(unique(bin_label))
        bin_label[findall(x -> x == el, bin_label)] .= idx
    end
    return res.binary_labels .= bin_label
end

"""
$(TYPEDSIGNATURES)

Removes unphysical solutions from the bitstrings in `res`.
"""
function clean_bitstrings(res::Result)
    return [
        [el for el in bit_string[phys_string]] for
        (bit_string, phys_string) in zip(res.classes["stable"], res.classes["physical"])
    ] #remove unphysical solutions from strings
end; #remove unphysical solutions from strings

"""
$(TYPEDSIGNATURES)

Converts a bit array to an integer.
"""
function bitarr_to_int(arr)
    return sum(arr .* (2 .^ collect((length(arr) - 1):-1:0)))
end

"""
$(TYPEDSIGNATURES)

Removes all solution branches from `res` where NONE of the solution falls into `class`.
Typically used to filter out unphysical solutions to prevent huge file sizes.
"""
function filter_result!(res::Result, class::String)
    bools = [any(getindex.(res.classes[class], i)) for i in 1:branch_count(res)]
    res.solutions = [s[bools] for s in res.solutions]
    for c in keys(res.classes)
        res.classes[c] = [s[bools] for s in res.classes[c]]
    end
end

"""
$(TYPEDSIGNATURES)

Classifies solutions in `result` using default classes: "physical", "stable, "Hopf".
"""
function _classify_default!(result)
    classify_solutions!(result, _is_physical, "physical")
    classify_solutions!(result, _is_stable(result), "stable")
    classify_solutions!(result, _is_Hopf_unstable(result), "Hopf")
    order_branches!(result, ["physical", "stable"]) # shuffle the branches to have relevant ones first
    classify_binaries!(result) # assign binaries to solutions depending on which branches are stable
    return nothing
end
