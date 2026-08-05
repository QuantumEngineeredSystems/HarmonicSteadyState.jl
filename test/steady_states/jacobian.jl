using HarmonicBalance, HarmonicSteadyState, OrderedCollections, Test
using Symbolics: Num
HB = HarmonicBalance
HSS = HarmonicSteadyState

# Stability is read off the implicit Jacobian. It must agree with the Jacobian of the
# rearranged system at every steady state, which is the only place it is evaluated: away
# from one the two differ by ∂(M⁻¹)/∂u ⋅ f(u), which vanishes when f(u) = 0.
#
# Only systems whose rearranged Jacobian is small enough to compile belong here. The one in
# issue #20 needs ten minutes of LLVM for that matrix, which is the reason the solver no
# longer builds it; `solves without the rearranged Jacobian` below covers that case instead.
function worst_disagreement(heq, swept, fixed)
    sw, fx = HSS.promote_types(OrderedDict(swept), OrderedDict(fixed))
    vars = HSS._free_symbols(heq, sw)
    rearranged = HSS.compile_matrix(HB.get_Jacobian(heq), vars; rules=fx)
    implicit = HSS.get_implicit_Jacobian(heq; sym_order=vars, rules=fx)

    res = get_steady_states(heq, sw, fx; show_progress=false)
    keys_ = collect(keys(sw))
    worst = 0.0
    for (i, point) in enumerate(res.solutions), solution in point
        all(isfinite, solution) || continue # solution arrays are NaN-padded
        all(isapprox.(imag.(solution), 0; atol=1e-8)) || continue
        v = ComplexF64[solution..., (sw[k][i] for k in keys_)...]
        worst = max(worst, maximum(abs.(rearranged(v) .- implicit(v))))
    end
    return worst
end

@testset "implicit Jacobian matches the rearranged one" begin
    @variables α ω ω0 F γ η λ ψ t x(t)

    duffing = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
    )
    add_harmonic!(duffing, x, ω)
    @test worst_disagreement(
        get_harmonic_equations(duffing),
        ω => collect(range(0.9, 1.1, 20)),
        (ω0 => 1.0, α => 1.0, γ => 0.01, F => 0.01),
    ) < 1e-12

    parametron = DifferentialEquation(
        d(x, t, 2) + (ω0^2 - λ * cos(2 * ω * t)) * x + α * x^3 + γ * d(x, t) ~
            F * cos(ω * t + ψ),
        x,
    )
    add_harmonic!(parametron, x, ω)
    @test worst_disagreement(
        get_harmonic_equations(parametron),
        ω => collect(range(0.9, 1.1, 20)),
        (ω0 => 1.0, α => 1.0, γ => 0.01, λ => 0.05, F => 0.01, ψ => 0.0),
    ) < 1e-12

    # Nonlinear damping, the term that blows the rearrangement up once several harmonics
    # are involved. At one harmonic it still compiles, so the agreement is checkable.
    nonlinear_damping = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) + η * x^2 * d(x, t) ~ F * cos(ω * t), x
    )
    add_harmonic!(nonlinear_damping, x, ω)
    @test worst_disagreement(
        get_harmonic_equations(nonlinear_damping),
        ω => collect(range(0.9, 1.1, 20)),
        (ω0 => 1.0, α => 1.0, γ => 0.01, η => 0.01, F => 0.01),
    ) < 1e-12
end

@testset "padded solutions give a NaN Jacobian, not an error" begin
    # `pad_solutions` fills short parameter points with NaN, and `_is_stable` evaluates the
    # Jacobian at those entries. The factorisation must return NaN there rather than throw.
    @variables α ω ω0 F γ t x(t)
    duffing = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
    )
    add_harmonic!(duffing, x, ω)
    heq = get_harmonic_equations(duffing)
    sw, fx = HSS.promote_types(
        OrderedDict(ω => collect(range(0.9, 1.1, 5))),
        OrderedDict(ω0 => 1.0, α => 1.0, γ => 0.01, F => 0.01),
    )
    jac = HSS.get_implicit_Jacobian(heq; sym_order=HSS._free_symbols(heq, sw), rules=fx)
    padded = ComplexF64[NaN, NaN, 1.0]
    @test all(isnan, jac(padded))
end

@testset "solves without the rearranged Jacobian" begin
    # Issue #20: nonlinear damping over two drive tones. The rearranged Jacobian is 366698
    # expression nodes and needs ~10 minutes to compile on first call; the implicit one is
    # 1036 nodes. This must stay well under a second once the machinery is warm.
    @variables α Δ ω0 δ F γ η t x(t)
    diff_eq = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) + η * x^2 * d(x, t) ~
            F * (cos((ω0 + δ - Δ / 2) * t) + cos((ω0 + δ + Δ / 2) * t)),
        x,
    )
    add_harmonic!(diff_eq, x, ω0 + δ - Δ / 2)
    add_harmonic!(diff_eq, x, ω0 + δ + Δ / 2)

    harmonic_eq = get_harmonic_equations(diff_eq)
    swept = δ => collect(range(-0.0625, 0.0125, 10))
    fixed = (Δ => 1.25, γ => 0.002, ω0 => 41.140, α => -3.2, η => 0.0078, F => 20.0)

    result = get_steady_states(harmonic_eq, swept, fixed; show_progress=false)
    @test length(result.solutions[1]) == 9
    @test sum(any.(get_class(result, "stable"))) == 2

    elapsed = @elapsed get_steady_states(harmonic_eq, swept, fixed; show_progress=false)
    @test elapsed < 10
end

@testset "autodiff backend gives the same Jacobian" begin
    using DifferentiationInterface: AutoForwardDiff
    import ForwardDiff

    @variables α ω ω0 F γ η t x(t)
    duffing = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) + η * x^2 * d(x, t) ~ F * cos(ω * t), x
    )
    add_harmonic!(duffing, x, ω)
    heq = get_harmonic_equations(duffing)
    swept = ω => collect(range(0.9, 1.1, 20))
    fixed = (ω0 => 1.0, α => 1.0, γ => 0.01, η => 0.01, F => 0.01)

    sw, fx = HSS.promote_types(OrderedDict(swept), OrderedDict(fixed))
    vars = HSS._free_symbols(heq, sw)
    implicit = HSS.get_implicit_Jacobian(heq; sym_order=vars, rules=fx)
    ad = HSS.get_ad_Jacobian(heq; sym_order=vars, rules=fx, backend=AutoForwardDiff())

    res = get_steady_states(heq, sw, fx; show_progress=false)
    worst = 0.0
    for (i, point) in enumerate(res.solutions), solution in point
        all(isfinite, solution) || continue # solution arrays are NaN-padded
        v = ComplexF64[solution..., sw[ω][i]]
        worst = max(worst, maximum(abs.(implicit(v) .- ad(v))))
    end
    @test worst < 1e-12 # holds at complex solutions too, the equations being holomorphic
    @test all(isnan, ad(ComplexF64[NaN, NaN, 1.0]))

    # `classify_solutions` evaluates the Jacobian under `Threads.@threads`, so a backend
    # holding shared mutable state misclassifies a solution now and then rather than always.
    reference = get_class(res, "stable")
    for _ in 1:5
        r = get_steady_states(
            heq, sw, fx; show_progress=false, jacobian_backend=AutoForwardDiff()
        )
        @test get_class(r, "stable") == reference
    end
end
