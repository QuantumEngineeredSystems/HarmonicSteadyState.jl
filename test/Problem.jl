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

    prob = HarmonicSteadyState.HomotopyContinuationProblem(
        [x, y], Num[], OrderedDict{Num,Vector{Float64}}(), OrderedDict{Num,Float64}(), F
    )
    @test_throws UndefRefError prob.jacobian
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
