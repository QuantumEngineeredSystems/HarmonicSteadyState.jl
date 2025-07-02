"""
Make a matrix that transforms the Jacobian to the dynamical matrix.

# Arguments
- `N::Int`: Size of system

# Returns
- `S, Sinv`: transformation matrices
"""
function make_S(N::Int)
    @assert iseven(N)
    k = N ÷ 2
    S_block = [1 im; 1 -im] / √(2)
    Sinv_block = [1 1; -im im] / √(2)

    S = kron(I(k), S_block)
    Sinv = kron(I(k), Sinv_block)
    return S, Sinv
end

"""
NO LONGER USED
Construct compiled function for the dynamical matrix of the system.

# Arguments
- `eqs::MeanfieldEquations`: Mean field equations of the system
- `problem::HarmonicSteadyState.HomotopyContinuationProblem`: HomotopyContinuationProblem of the system

# Returns
- `Dnum`: Compiled function to calculate dynamical matrix of the system

"""
function make_D(eqs::HarmonicEquation, varied, fixed)
    jac = eqs.jacobian

    variables = QuestBase._remove_brackets(QuestBase.get_variables(eqs))
    S, Sinv = make_S(length(variables))
    Duv = expand.(S * jac * Sinv)

    subs = Num[variables..., keys(Dict(varied))...]
    Dtemp = Symbolics.substitute.(Duv, Ref(Dict(fixed)))
    jacfunc = Symbolics.build_function(Dtemp, subs; expression=Val(false))
    jacfunc = jacfunc isa Tuple ? jacfunc[1] : jacfunc

    return HarmonicSteadyState.JacobianFunction(ComplexF64)(jacfunc)
end

"""
$(TYPEDSIGNATURES)

Compute the forward_transmission response spectrum, i.e, how much of an input signal applied
at port 1 emerges at port 2. Colloquially known as S21 parameter in microwave engineering.
The amplitude and phase of S21 tell you how much signal is transmitted and with what delay
or phase shift:
- If S21 ≈ 1 (or 0 dB), the system transmits all power from input to output.
- If S21 ≈ 0 (or very negative dB), very little signal is transmitted.

# Arguments
- `result::Result`: Result object containing the system's solutions
- `op_index::Int=1`: Index of operator in mean field equations to evaluate response for
- `Ω_range`: Range of frequencies to evaluate
- `branch::Int`: Branch number to analyze
- `class="stable"`: Class of solutions to evaluate response for

# Returns
- `χ`: Complex response matrix where rows correspond to frequencies and columns to solutions

# Example
```julia
Ω_range = range(-0.2, 0.2, 500)

χ = get_susceptibility(result, Ω_range, 3 #=branch=#);

κ_ext = 0.05
S21 = 1 .- χ*κ_ext/2
S21_log = 20 .* log10.(abs.(S21)) # expressed in dB
```
"""
function get_susceptibility(
    result::HarmonicSteadyState.Result, op_index::Int, Ω_range, branch::Int; class="stable"
)
    vars = result.problem.variables
    S, Sinv = make_S(length(vars))

    inclass = get_class(result, branch, class)
    !any(inclass) && error("Cannot generate a spectrum - no $class solutions!")

    inds = findall(inclass)

    χ = Matrix{ComplexF64}(undef, length(Ω_range), length(inds))

    Symbolics.@variables Ω

    for (i, index) in enumerate(inds)
        sol = get_single_solution(result; branch, index)

        jac = result.jacobian(collect(values(sol)))

        D = S * jac * Sinv

        eig = eigen(D)

        χ_mat =
            -eig.vectors *
            LinearAlgebra.Diagonal(1 ./ (eig.values .- im * Ω)) *
            LinearAlgebra.inv(eig.vectors)

        χ_func = Symbolics.build_function(
            χ_mat[op_index, op_index], Ω; expression=Val{false}
        )

        for (j, Ω) in enumerate(Ω_range)
            @inbounds χ[j, i] = χ_func(Ω)
        end
    end
    return χ
end
