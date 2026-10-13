using HarmonicSteadyState
using HarmonicSteadyState: OrderedDict
using Test

@testset "Problem without EOM" begin
    using HomotopyContinuation: System, Expression
    using Symbolics: Num, @variables
    F = System([Expression("x^2")])

    @variables x y
    prob = HarmonicSteadyState.HomotopyContinuationProblem(
        [x, y],
        Num[],
        OrderedDict{Num,Vector{Float64}}(),
        OrderedDict{Num,Float64}(),
        F,
        HarmonicSteadyState.JacobianFunction(ComplexF64)(x -> x),
    )
    @test isnothing(HarmonicSteadyState.source(prob))
    @test HarmonicSteadyState.source_type(prob) == Nothing
    @test_throws ErrorException HarmonicSteadyState._harmonic_equation(prob, "do a thing")

    prob = HarmonicSteadyState.HomotopyContinuationProblem(
        [x, y], Num[], OrderedDict{Num,Vector{Float64}}(), OrderedDict{Num,Float64}(), F
    )
    @test_throws UndefRefError prob.jacobian
    @test isnothing(HarmonicSteadyState.source(prob))
    @test HarmonicSteadyState.source_type(prob) == Nothing
end

@testset "explicit symbol order with Unicode subscripts" begin
    using Symbolics: @variables
    @variables a₂ᵣ b₁ᵢ p₁ᵣ
    problem = HarmonicSteadyState.HomotopyContinuationProblem(
        [a₂ᵣ - b₁ᵢ^2, b₁ᵢ - p₁ᵣ],
        [a₂ᵣ, b₁ᵢ],
        [p₁ᵣ],
        OrderedDict(p₁ᵣ => [1.0, 2.0]),
        OrderedDict{Num,Float64}(),
    )
    @test string.(problem.variables) == ["a₂ᵣ", "b₁ᵢ"]
    @test string.(problem.parameters) == ["p₁ᵣ"]
    result = get_steady_states(
        problem, TotalDegree(); sorting=:none, classify_default=false, show_progress=false
    )
    for (solutions, expected) in zip(result.solutions, ([1.0, 1.0], [4.0, 2.0]))
        finite = filter(solution -> all(isfinite, solution), solutions)
        @test !isempty(finite)
        @test all(solution -> isapprox(solution, expected; atol=1e-6), finite)
    end
end

@testset "Problem with self made system" begin
    using Symbolics
    @variables u1 v1 Δ F
    eqs = [Δ * u1 + F, -Δ * v1]
    vars = [u1, v1]
    pars = [Δ, F]
    fix = OrderedDict{Num,Float64}(F => 0.005)
    swept = OrderedDict{Num,Vector{Float64}}(Δ => range(0.0, 1.0; length=10))
    prob = HarmonicSteadyState.HomotopyContinuationProblem(eqs, vars, pars, swept, fix)
    @test all(string.(prob.variables) .== string.(vars))
    @test all(string.(prob.parameters) .== string.(pars))
    @test all(string.(prob.system.expressions) .== string.(eqs))
end
