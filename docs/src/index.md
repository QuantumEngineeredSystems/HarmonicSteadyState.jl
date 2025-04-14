# HarmonicSteadyState.j

A package for computing the classical steady state of the effective stroboscopic dynamical systems. Given one has the autonomous equations of motion of the system in the rotating frame of the characteristic response frequencies, it collect steady states methods to find and describe the stationary responses of the system. It supports the following methods:

- fixed point steady states with Homotopy Continuation
- Finding Limit-cycle with Homotopy Continuation
- Stability analysis
- Linear response of the steady state in the (non-)rotating frame
- Parameter sweeps
- Plotting utilities

## Other packages in the Quest Ecosystem compatible with HarmonicSteadyState.jl

### [HarmonicBalance.jl](https://github.com/QuantumEngineeredSystems/HarmonicBalance.jl)

A package for applying the harmonic balance method to classical driven nonlinear dynamical systems.
It computes the stroboscopic effective equations of motion of the system at the characteristic response frequencies of the system.
Both [Krylov-Bogoliubov](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/stable/manual/extracting_harmonics#Krylov-Bogoliubov) averaging method to higher orders and the [harmonic balance method](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/stable/manual/extracting_harmonics#Krylov-Bogoliubov) are implemented.
