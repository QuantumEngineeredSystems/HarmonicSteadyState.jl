# HarmonicSteadyState.jl

[![docsdev](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://quantumengineeredsystems.github.io/HarmonicBalance/dev
)
[![docsstable](https://img.shields.io/badge/docs-stable-blue.svg)](https://quantumengineeredsystems.github.io/HarmonicBalance/stable
)
[![codecov](https://codecov.io/gh/QuantumEngineeredSystems/HarmonicSteadyState/branch/main/graph/badge.svg)](https://codecov.io/gh/QuantumEngineeredSystems/HarmonicSteadyState.jl)

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%E2%9C%88%EF%B8%8F%20tested%20with%20-%20JET.jl%20-%20red)](https://github.com/aviatesk/JET.jl)

A package for computing the classical steady state of the effective stroboscopic dynamical systems. Given one has the autonomous equations of motion of the system in the rotating frame of the characteristic response frequencies, it collect steady states methods to find and describe the stationary responses of the system. It supports the following methods:

- fixed point steady states with Homotopy Continuation
- Finding Limit-cycle with Homotopy Continuation
- Stability analysis
- Linear response of the steady state in the (non-)rotating frame
- Parameter sweeps
- Plotting utilities
