# multiply a peak by a number.
#! format: off
function Base.:*(number::T, peak::Lorentzian{T}) where {T<:Real} # multiplication operation
    return Lorentzian(peak.ω0, peak.Γ, *(peak.A, number))
end
function Base.:*(peak::Lorentzian{T}, number::T) where {T<:Real} # multiplication operation
    return Lorentzian(peak.ω0, peak.Γ, *(peak.A, number))
end
#! format: on

# TODO:  ∨ allocation heavy. Instead, define in-place multipliation.
function Base.:*(number::T, s::JacobianSpectrum{T}) where {T<:Real}
    return JacobianSpectrum{T}([*(number, peak) for peak in s.peaks])
end
function Base.:*(s::JacobianSpectrum{T}, number::T) where {T<:Real}
    return JacobianSpectrum{T}([*(number, peak) for peak in s.peaks])
end

function Base.show(io::IO, s::JacobianSpectrum)
    peaks = sort(s.peaks; by=x -> x.ω0)
    print(io, "Lorentzian peaks (central frequency ω0, linewidth Γ): \n")
    for peak in peaks
        @printf(io, "%.3e * L(ω0 = %.6e, Γ = %.3e)\n", peak.A, peak.ω0, peak.Γ)
    end
end

function Base.show(io::IO, spectra::Dict{Num,JacobianSpectrum})
    for var in keys(spectra)
        print("VARIABLE: ", var, "\n")
        show(spectra[var])
        print("\n")
    end
end

"Adds a peak p to JacobianSpectrum s."
function add_peak(s::JacobianSpectrum, p::Lorentzian)
    return JacobianSpectrum(vcat(s.peaks, p))
end

"Adds all peaks from s2 to s1."
function add_peak(s1::JacobianSpectrum, s2::JacobianSpectrum)
    return JacobianSpectrum(vcat(s1.peaks, s2.peaks))
end

"Adds a peak `p` to `s`, growing `s` in place."
function add_peak!(s::JacobianSpectrum, p::Lorentzian)
    push!(s.peaks, p)
    return s
end

"Adds all peaks from `s2` to `s1`, growing `s1` in place."
function add_peak!(s1::JacobianSpectrum, s2::JacobianSpectrum)
    append!(s1.peaks, s2.peaks)
    return s1
end

"Gives the numerical value of a peak at ω."
function evaluate(peak::Lorentzian{T}, ω::T) where {T<:Real}
    return peak.A / sqrt(((peak.ω0 - ω)^2 + (peak.Γ)^2))
end

"Gives the numerical value of a JacobianSpectrum at ω"
function evaluate(s::JacobianSpectrum{T}, ω::T) where {T<:Real}
    sum = zero(T)
    for p in s.peaks
        sum += evaluate(p, ω)
    end
    return sum
end

"""
Take a pair of harmonic variable u,v and an eigenvalue λ and eigenvector eigvec_2d
of the Jacobian to generate corresponding Lorentzians.
IMPORTANT: The eigenvetor eigen_2d contains only the two components of the full eigenvector
which correspond to the u,v pair in question.
"""
function _pair_to_peaks(λ, eigvec_2d::Vector; ω)
    u, v = eigvec_2d
    peak1 =
        (1 + 2 * (imag(u) * real(v) - real(u) * imag(v))) *
        Lorentzian(; ω0=ω - imag(λ), Γ=real(λ))
    peak2 =
        (1 + 2 * (real(u) * imag(v) - real(v) * imag(u))) *
        Lorentzian(; ω0=ω + imag(λ), Γ=real(λ))
    return JacobianSpectrum([peak1, peak2])
end

"Return the indices of a-type variables (zero harmonics) from hvars."
_get_as(hvars::Vector{HarmonicVariable}) = findall(x -> isequal(x.type, "a"), hvars)

#   Returns the spectra of all variables in `res` for `index` of `branch`.
"""
Here linear response is treated with the slow-flow approximation (SFA), see Chapter 5
of JK's thesis. Linear response always appears as a sum of Lorentzians, but is inaccurate where
these are peaked far from the drive frequency.
"""
function JacobianSpectrum(
    res::Result{D,S,P}; index::Int, branch::Int, force=false
) where {D,S,P}
    eom = QuestBase.source(res.problem)
    hvars = eom.variables # fetch the vector of HarmonicVariable
    # blank JacobianSpectrum for each variable
    all_spectra = Dict{Num,JacobianSpectrum{P}}(
        hvar.natural_variable => JacobianSpectrum{P}() for hvar in hvars
    )

    if force
        res.classes["stable"][index][branch] || return all_spectra # if the solution is unstable, return empty spectra
    else
        res.classes["stable"][index][branch] ||
            error("\nThe solution is unstable - it has no JacobianSpectrum!\n")
    end

    solution_dict = get_single_solution(res; branch=branch, index=index)
    solutions = get_variable_solutions(res; branch=branch, index=index)
    λs, vs = eigen(res.jacobian(solutions))

    uv_pairs = _get_uv_pairs(hvars)
    a_idxs = _get_as(hvars)

    # Neither the harmonic of a uv pair nor the spectrum a variable accumulates into depends
    # on the eigenvalue, so both are resolved once per variable instead of once per
    # (eigenvalue, variable). The symbolic substitution dominates this loop, so hoisting it
    # matters.
    uv_ωnums = [
        real(
            SymbolicUtils.unwrap_const(
                SymbolicUtils.unwrap(
                    Symbolics.substitute(hvars[first(pair)].ω, solution_dict)
                ),
            ),
        ) for pair in uv_pairs
    ]
    uv_spectra = [all_spectra[hvars[first(pair)].natural_variable] for pair in uv_pairs]
    a_spectra = [all_spectra[hvars[a_idx].natural_variable] for a_idx in a_idxs]

    for (j, λ) in enumerate(λs)
        eigvec = @view vs[:, j] # the eigenvector

        # 2 peaks for each pair of uv variables
        for (p, pair) in enumerate(uv_pairs)
            eigvec_2d = eigvec[pair] # fetch the relevant part of the Jacobian eigenvector
            ωnum = uv_ωnums[p]
            # ^ the harmonic (numerical now) associated to this harmonic variable

            # eigvec_2d is associated to a natural variable -> this variable gets Lorentzian peaks
            peaks = norm(eigvec_2d) * _pair_to_peaks(λ, eigvec_2d; ω=ωnum)

            add_peak!(uv_spectra[p], peaks)
        end

        # 1 peak for a-type variable
        for (p, a_idx) in enumerate(a_idxs)
            eigvec_1d = eigvec[a_idx]
            peak = 2 * norm(eigvec_1d) * Lorentzian(; ω0=abs(imag(λ)), Γ=real(λ))
            add_peak!(a_spectra[p], peak)
        end
    end
    #_simplify_spectra!(all_spectra) # condense similar peaks
    return all_spectra
end

"Condense peaks with similar centers and linewidths in JacobianSpectrum s."
function _simplify_spectrum(s::JacobianSpectrum{T}) where {T<:Real}
    ω0s = getfield.(s.peaks, Symbol("ω0"))
    new_spectrum = JacobianSpectrum{T}()
    for ω in unique(ω0s)
        peaks = filter(x -> x.ω0 == ω, s.peaks) # all peaks centered at ω
        length(unique(getfield.(peaks, Symbol("Γ")))) > 1 && return (s)  # there are peaks of multiple linewidths at the same frequency => not simplifiable
        total_height = sqrt(sum(getfield.(peaks, Symbol("A")) .^ 2)) # new height is the sum of squares
        new_spectrum = add_peak(
            new_spectrum, total_height * Lorentzian(; ω0=ω, Γ=peaks[1].Γ)
        )
    end
    return new_spectrum
end

"Simplify every JacobianSpectrum in a dictionary assigning variables to spectra."
function _simplify_spectra!(spectra::Dict{Num,JacobianSpectrum})
    for var in keys(spectra)
        spectra[var] = _simplify_spectrum(spectra[var])
    end
end
