# New functions

using LinearAlgebra
using HarmonicSteadyState.Symbolics
QCExt = Base.get_extension(HarmonicSteadyState, :QuantumCumulantsExt)

"""
Make a matrix that transforms the Jacobian to the dynamical matrix.

# Arguments
- `N::Int`: Size of system

# Returns
- `S, Sinv`: transformation matrices
"""
function make_S_Sinv(N::Int)
	@assert iseven(N)
	k = N ÷ 2
	S_block = [1 im; 1 -im] / √(2)
	Sinv_block = [1 1; -im im] / √(2)

	S = kron(I(k), S_block)
	Sinv = kron(I(k), Sinv_block)
	return S, Sinv
end

"""
Construct compiled function for the dynamical matrix of the system.

# Arguments
- `eqs::MeanfieldEquations`: Mean field equations of the system
- `problem::HarmonicSteadyState.HomotopyContinuationProblem`: HomotopyContinuationProblem of the system

# Returns
- `Dnum`: Compiled function to calculate dynamical matrix of the system
"""
function make_D(eqs::MeanfieldEquations, problem::HarmonicSteadyState.HomotopyContinuationProblem)
	variables, equations = QCExt.compute_real_equations(eqs)
	jac = HarmonicSteadyState.get_Jacobian(equations, variables)

	S, Sinv = make_S_Sinv(length(variables))

	Duv = S * jac * Sinv .|> expand

	varied = problem.swept_parameters
	fixed = problem.fixed_parameters

	subs = Num[variables..., keys(Dict(varied))...]
	Jtemp = Symbolics.substitute.(Duv, Ref(Dict(fixed)))

	jacfunc = Symbolics.build_function(Jtemp, subs; expression = Val(false))

	Dnum = HarmonicSteadyState.JacobianFunction(ComplexF64)(jacfunc isa Tuple ? jacfunc[1] : jacfunc)
	return Dnum
end

"""
Calculate the Jacobian response spectrum for a system construced from a Hamiltonian using Quantum get_cumulant_response.

# Arguments
- `eqs::MeanfieldEquations`: Mean field equations of the system
- `res::Result`: Result object containing the system's solutions
- `Ω_range`: Range of frequencies to evaluate
- `branch::Int`: Branch number to analyze
- `op_index::Int=1`: Index of operator in mean field equations to evaluate response for
- `use_stable=true`: Evaluate response only for stable steady states, or also for any physical solutions

# Returns
- `χ`: Complex response matrix where rows correspond to frequencies and columns to solutions
"""
function get_linear_response_cumulants(eqs::MeanfieldEquations, res::HarmonicSteadyState.Result, Ω_range, branch::Int; op_index::Int = 1, use_stable = true)
	D_func = make_D(eqs, res.problem)

	if use_stable
		stable = get_class(res, branch, "stable")
		!any(stable) && error("Cannot generate a spectrum - no stable solutions!")
	else
		stable = get_class(res, branch, "physical")
	end

	χ = Matrix{ComplexF64}(undef, length(Ω_range), length(findall(stable)))

	@variables Ω

	for (i, ind) in enumerate(findall(stable))
		sol = HarmonicSteadyState.get_variable_solutions(res; branch = branch, index = ind)

		D = D_func(real.(sol))
		eig = eigen(D)

		χ_mat = -eig.vectors * Diagonal(1 ./ (eig.values .- im*Ω)) * inv(eig.vectors)

		χ_func = Symbolics.build_function(χ_mat[op_index, op_index], Ω; expression = Val{false})

		for (j, Ω) in enumerate(Ω_range)
			@inbounds χ[j, i] = χ_func(Ω)
		end
	end
	return χ
end
