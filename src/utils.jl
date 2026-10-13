function common_type(A::Array)
    types = unique(reduce(vcat, [unique(typeof.(ps)) for ps in A]))
    return promote_type(types...)
end
function collect_parameters(sweeps::OrderedDict, fixed_parameters::OrderedDict)
    param_ranges = collect(values(sweeps))
    iter = Iterators.product(param_ranges..., values(fixed_parameters)...)
    return collect(iter)
end
function type_stable_parameters(sweeps::OrderedDict, fixed_parameters::OrderedDict)
    input_array = collect_parameters(sweeps, fixed_parameters)
    T = common_type(input_array)
    return [convert.(T, ps) for ps in input_array]
end
function parameter_type(sweeps::AbstractDict, fixed_parameters::AbstractDict)
    T = Union{}
    for range in values(sweeps), value in range
        T = promote_type(T, typeof(value))
    end
    for value in values(fixed_parameters)
        T = promote_type(T, typeof(value))
    end
    return T
end

function promote_types(sweeps::AbstractDict, fixed_parameters::AbstractDict)
    T = parameter_type(sweeps, fixed_parameters)

    promoted_sweeps = OrderedDict{Num,Vector{T}}(
        key => Vector{T}(range) for (key, range) in pairs(sweeps)
    )
    promoted_fixed = OrderedDict{Num,T}(fixed_parameters)
    return promoted_sweeps, promoted_fixed
end

"Remove occurrences of `sweeps` elements from `fixed_parameters`."
function filter_duplicate_parameters(sweeps, fixed_parameters)
    new_params = copy(fixed_parameters)
    for par in keys(sweeps)
        delete!(new_params, par)
    end
    return new_params
end

""" Project the array `a` into the real axis, warning if its contents are complex. """
function _realify(a::Array{T}; warning="") where {T<:Number}
    warned = false
    a_real = similar(a, typeof(real(a[1])))
    for i in eachindex(a)
        if !isnan(a[i]) && !warned && !is_real(a[i])
            @warn "Values with non-negligible complex parts have
            been projected on the real axis! " * warning
            warned = true
        end
        a_real[i] = real(a[i])
    end
    return a_real
end

_realify(a::Array{Real}; warning="") = a

is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL || abs(x) < 1e-70
is_real(x::Array) = is_real.(x)
is_zero(x::Vector{ComplexF64}; atol=1e-6) = all(isapprox.(x, complex(0.0); atol))

_str_to_vec(s::Vector) = s
_str_to_vec(s) = [s]

_parse_expression(x::Expr)::Num = Symbolics.parse_expr_to_symbolic(x, @__MODULE__)
_parse_expression(x::Symbol)::Num = Symbolics.parse_expr_to_symbolic(x, @__MODULE__)
_parse_expression(x::String)::Num = _parse_expression(Meta.parse(x))

tuple_to_vector(t::Tuple) = [i for i in t]
