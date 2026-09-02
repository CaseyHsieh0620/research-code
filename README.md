# research

A few pieces of code from my open quantum systems work.

## redfield-solver

A C++ implementation of the Redfield equation, meant to simulate quantum dynamics with non-Markovian memory. It builds the system Hamiltonian and system-bath coupling, rotates into the energy eigenbasis, computes the Redfield dissipation tensor from an Ohmic spectral density, and propagates the reduced density matrix forward in time.

Inputs are plain text: `Hamiltonian.txt` for the system Hamiltonian, `Sm.txt` for the system-bath coupling operator, and `parameters.txt` for the bath and propagation settings (temperature, cutoff frequency, coupling strength, timestep, number of steps). `parameters.cpp` reads these into a `BathParam` struct that the solver runs on. Output is the full time series of reduced density matrix elements.

## level-statistics

A generator of randomized Hamiltonians, used to study eigenvalue statistics and level spacings.

## coherence-map-visualizer

A tool built for the Makri group, meant to integrate with her PATHSUM code to build coherence maps. It takes the reduced density matrix output from a PATHSUM run, pulls out the matrix at whatever time you ask for, and plots the real and imaginary parts side by side so you can actually see the coherences instead of squinting at a column of numbers.
