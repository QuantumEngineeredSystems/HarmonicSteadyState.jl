```@meta
EditURL = "../../../examples/forward_transmission.jl"
```

````@example forward_transmission
# Linear response and transmission/reflection coeffictients for magnon three-wave mixing

using HarmonicSteadyState, QuantumCumulants, Plots
````

Consider a model of nonlinear magnon-magnon coupling, as described in [this](https://arxiv.org/abs/2506.11527) paper. The model describes  a three-wave mixing interaction between a strongly driven $k=0$ FMR mode and two parametricallly excited propagating modes with opposite momentum $\pm k$ and at half frequency. In this notebook, we will show how the cumulant approximation can be used to calculate steady states, and to calculate the $S_{21}$ transmission coefficient.

````@example forward_transmission
hm = FockSpace(:magnon)
hc = FockSpace(:polariton)
h = hm ⊗ hc  # Hilbertspace

@qnumbers m::Destroy(h, 1) c::Destroy(h, 2) # Operators

@rnumbers Δ Vk Ωd γm γk # Parameters
param = [Δ, Vk, Ωd, γm, γk]

H_RWA_sym = (
    Δ * m' * m + Δ / 2 * c' * c + Vk * m * c' * c' + Vk * m' * c * c + (Ωd * m + Ωd * m')
)
ops = [m, m', c, c'] # Operators for meanfield evolution

eqs_RWA = meanfield(ops, H_RWA_sym, [m, c]; rates=[γm, γk], order=1)
eqs_completed_RWA = complete(eqs_RWA) # Meanfield equations
````

Let's sweep the power of the drive $\Omega_d$, with $\Delta=0$, and solve for the steady state. The steady-state solutions show that the FMR mode saturates after a threshold power, followed by the coherent excitation of the parametrically induced counter-propagating modes.

````@example forward_transmission
drive_range = range(0, 1, 100)

fixed = (Δ => 0, Vk => 0.0002, γm => 0.1, γk => 0.01)
varied = (Ωd => drive_range)
problem = HarmonicSteadyState.HomotopyContinuationProblem(
    eqs_completed_RWA, param, varied, fixed
)

result = get_steady_states(problem, TotalDegree())

plot(plot(result; y="1/sqrt(2)*(mᵣ+ mᵢ)"), plot(result; y="1/sqrt(2)*(cᵣ + cᵢ)"))

# Linear response and S21
````

To find the response of the driven system to a second, weak probe, we use the method described [here](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/stable/background/stability_response#linresp_background). Here, we calculate the response in the same rotating frame as the Hamiltonian. The linear response is related to the scattering parameter $S_{21}$ by
$$S_{21}=1-\frac{\kappa_{ext}}{2} \chi,$$
where $\kappa_{ext}$ is the coupling of the system to the measurement apparatus.

The result below shows the characteristic splitting of the magnon resonance above the power threshold, which matches the experiment.

````@example forward_transmission
Ω_range = range(-0.2, 0.2, 500)
χ = get_forward_transmission_response(result, Ω_range, 3);

κ_ext = 0.05
S21 = 1 .- χ*κ_ext/2
S21_log = 20 .* log10.(abs.(S21))

stable = get_class(result, 3, "physical")
y = drive_range[stable]
heatmap(Ω_range, y, transpose(S21_log); c=:matter, cbar_title="S21 (dB)")
ylabel!("Ω_d")
xlabel!("Probe detuning")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

