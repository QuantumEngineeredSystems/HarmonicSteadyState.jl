"""
Compile the Jacobian from `prob`, inserting `fixed_parameters`.
Returns a function that takes a vector of variables and `swept_parameters` to give the Jacobian.
The order of the vector is first the variables, then the swept parameters.
For LC variables called "hopf", the Jacobian is already compiled.
For Type stability the function is wrapped in a FunctionWrapper.

The Jacobian is evaluated implicitly, see [`get_implicit_Jacobian`](@ref), or with `backend`
if one is given, see [`get_ad_Jacobian`](@ref). Neither reads `eom.jacobian`, the Jacobian of
the system rearranged so the derivatives stand alone on one side. All three give the same
matrix at a steady state, but the rearrangement inverts the mass matrix symbolically and the
result is compiled as one expression:

| system                   | rearranged | implicit | first call | per call    |
|:-------------------------|-----------:|---------:|-----------:|------------:|
| Duffing                  |    181     |    69    | 0.02/0.13s | 0.17/1.18µs |
| two-tone Duffing         |   1632     |   340    | 0.10/0.13s | 0.60/2.38µs |
| two-tone, `η*x^2*d(x,t)` | 366698     |  1036    | 583/0.12s  | 610/5.9µs   |

Nonlinear damping over several harmonics blows the rearranged expression up by two orders of
magnitude, and LLVM then needs ten minutes on the first call. The implicit route pays a small
dense solve per call instead, which costs about a microsecond on the systems where the
rearranged Jacobian was cheap and wins outright everywhere else.
"""
function _compile_Jacobian(
    eom::HarmonicEquation,
    soltype::DataType,
    swept::OrderedDict,
    fixed::OrderedDict;
    backend=nothing,
)::JacobianFunction(soltype)
    sym_order = _free_symbols(eom, swept)
    compiled_J = if "Hopf" ∈ getfield.(eom.variables, :type)
        eom.jacobian # already a compiled gauge-fixed function, see LimitCycles
    elseif isnothing(backend)
        get_implicit_Jacobian(eom; sym_order, rules=fixed)
    else
        get_ad_Jacobian(eom; sym_order, rules=fixed, backend)
    end
    return JacobianFunction(soltype)(compiled_J)
end

"""
Take a matrix containing symbolic variables `variables` and keys of `fixed_parameters`.
Substitute the values according to `fixed_parameters` and compile into a function that takes
numerical arguments in the order set in `variables`.
"""
function compile_matrix(
    mat::Matrix{Num}, variables::Vector{Num}; rules=Dict()
)::RuntimeGeneratedFunction
    J = substitute_all.(mat, Ref(rules)) # Ref makes sure only mat is broadcasted
    # `cse` shares repeated subexpressions across the whole matrix, which is what the entries
    # of a Jacobian mostly are. It cuts the LLVM compile on the first call by 2-3x on larger
    # systems and the per-call cost by about a tenth, to the last bit of the same result.
    jacfunc = Symbolics.build_function(J, variables; expression=Val(false), cse=true)
    return jacfunc isa Tuple ? jacfunc[1] : jacfunc
end

function compute_and_compile_Jacobian(
    equations::Vector{Num},
    variables::Vector{Num},
    soltype::DataType,
    swept::AbstractDict,
    fixed::AbstractDict,
)::JacobianFunction(soltype)
    subs = Num[variables..., keys(swept)...]
    # jac = Symbolics.jacobian(equations, variables)
    jac = get_Jacobian(equations, variables)
    return JacobianFunction(soltype)(compile_matrix(jac, subs; rules=fixed))
end

function _free_symbols(eom::HarmonicEquation, swept::OrderedDict)::Vector{Num}
    return cat(declare_variables(eom), collect(keys(swept)); dims=1)
end

"""
Code follows for an implicit treatment of the Jacobian. Usually we rearrange the linear response
equations to have time-derivatives on one side. This may be extremely costly. Implicit evaluation
means only solving the equations AFTER numerical values have been plugged in, giving a constant
time cost per run.
"""

# # TODO COMPILE THIS?
"""
$(TYPEDSIGNATURES)

Construct a function for the Jacobian of `eom` using `rules=Dict()`.

Necessary matrix inversions are only performed AFTER substituting numerical values at each call,
avoiding huge symbolic operations.

Returns a function `f(soln::OrderedDict{Num,T})::Matrix{T}`.
"""
function get_implicit_Jacobian(eom::HarmonicEquation; sym_order, rules=Dict())
    J0c = compile_matrix(_get_J_matrix(eom; order=0), sym_order; rules)
    J1c = compile_matrix(_get_J_matrix(eom; order=1), sym_order; rules)
    function jacfunc(vals::Vector)
        # Both blocks are freshly allocated by their compiled functions, so factorising and
        # solving in place is safe and cuts the per-call cost sevenfold over `-inv(J1)*J0`.
        J1 = real.(J1c(vals))
        J0 = J0c(vals)
        # `check=false` to match `inv`/`\`: solution arrays are NaN-padded, and a Jacobian
        # at a padded entry has to come back NaN rather than throw.
        ldiv!(lu!(J1; check=false), J0)
        return J0 .= .-J0
    end
    return jacfunc
end

"Compile the harmonic equations as `f([u..., u̇..., swept...])`, the residual the Jacobian blocks are the derivatives of."
function _compile_residual(eom::HarmonicEquation, sym_order::Vector{Num}; rules=Dict())
    T = get_independent_variables(eom)[1]
    vars = get_variables(eom)
    n = length(vars)
    dsyms = [declare_variable("__du" * string(i)) for i in 1:n]
    subs = merge(Dict(d(vars, T, 1) .=> dsyms), Dict(vars .=> _remove_brackets.(vars)))
    R = Num[
        substitute_all(Symbolics.expand_derivatives(eq.lhs - eq.rhs), subs) for
        eq in eom.equations
    ]
    R = substitute_all.(R, Ref(rules))
    args = Num[sym_order[1:n]..., dsyms..., sym_order[(n + 1):end]...]
    f = Symbolics.build_function(R, args; expression=Val(false), cse=true)
    return (f isa Tuple ? f[1] : f), n
end

_split(o) = vcat(real.(o), imag.(o))

"""
$(TYPEDSIGNATURES)

Construct a function for the Jacobian of `eom`, differentiating the harmonic equations with
`backend` rather than symbolically.

`backend` is an ADTypes backend such as `AutoForwardDiff()`; load the differentiation package
it names to make it available. Both this and [`get_implicit_Jacobian`](@ref) solve for the
derivatives only once the parameters have values, so they agree at every steady state.

Returns a function `f(vals::Vector)::Matrix`.
"""
function get_ad_Jacobian(eom::HarmonicEquation; sym_order, rules=Dict(), backend)
    f, n = _compile_residual(eom, sym_order; rules)

    # The harmonic equations are polynomial, hence holomorphic, so a perturbation along the
    # real axis already carries the full complex derivative. Differentiating that and splitting
    # the result into real and imaginary parts keeps every backend on the real inputs they all
    # accept, while still admitting the complex solutions the Jacobian is evaluated at.
    gu(a, c) = _split(f(vcat(complex.(a, c[1]), zeros(Complex{eltype(a)}, n), c[2])))
    gd(b, c) = _split(f(vcat(c[1], complex.(b, zeros(length(b))), c[2])))

    # No `prepare_jacobian`: `classify_solutions` evaluates the Jacobian under `Threads.@threads`,
    # and a preparation object holds buffers the backend writes into, so sharing one across
    # threads races and silently misclassifies the odd solution.
    function jacfunc(vals::Vector)
        u = ComplexF64.(vals[1:n])
        rest = ComplexF64.(vals[(n + 1):end])
        M0 = DifferentiationInterface.jacobian(
            gu, backend, real.(u), Constant((imag.(u), rest))
        )
        M1 = DifferentiationInterface.jacobian(
            gd, backend, zeros(n), Constant((u, rest))
        )
        J0 = complex.(M0[1:n, :], M0[(n + 1):(2n), :])
        J1 = M1[1:n, :] # real part: the mass matrix has no imaginary component
        # `check=false` to match the implicit route: solution arrays are NaN-padded, and a
        # Jacobian at a padded entry has to come back NaN rather than throw.
        ldiv!(lu!(J1; check=false), J0)
        return J0 .= .-J0
    end
    return jacfunc
end

function get_implicit_Jacobian(p::SteadyStateProblem)
    eom = _harmonic_equation(p, "compile an implicit Jacobian")
    return get_implicit_Jacobian(eom; sym_order=_free_symbols(p), rules=p.fixed_parameters)
end

# for implicit evaluation, the numerical values precede the rearrangement
# for limit cycles, the zero eigenvalue causes the rearrangement to fail -> filter it out
# THIS SETS ALL DERIVATIVES TO ZERO - assumes use for steady states
function _get_J_matrix(eom::HarmonicEquation; order=0)
    order > 1 && error("Cannot get a J matrix of order > 1 from the harmonic equations.\n
                       These are by definition missing higher derivatives")

    vars_simp = Dict([var => _remove_brackets(var) for var in get_variables(eom)])
    T = get_independent_variables(eom)[1]
    # J = Symbolics.jacobian(eom.equations, d(get_variables(eom), T, order))
    J = get_Jacobian(eom.equations, d(get_variables(eom), T, order))

    return Symbolics.expand_derivatives.(substitute_all(J, vars_simp)) # a symbolic matrix to be compiled
end

# ∨ this is atm a temporary duplicate code from HarmonicBalance
" Get the Jacobian of a set of equations `eqs` with respect to the variables `vars`. "
function get_Jacobian(eqs::Vector{Num}, vars::Vector{Num})::Matrix{Num}
    length(eqs) == length(vars) || error("Jacobians are only defined for square systems!")
    M = Matrix{Num}(undef, length(vars), length(vars))

    for idx in CartesianIndices(M)
        M[idx] = Symbolics.expand_derivatives(d(eqs[idx[1]], vars[idx[2]]))
    end
    return M
end # TODO should replace with Symbolics.jacobian

function get_Jacobian(eqs::Vector{Symbolics.Equation}, vars::Vector{Num})::Matrix{Num}
    expr = Num[getfield(eq, :lhs) - getfield(eq, :rhs) for eq in eqs]
    return get_Jacobian(expr, vars)
end
