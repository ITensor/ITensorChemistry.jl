# ITensorChemistry

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mtfishman.github.io/ITensorChemistry.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mtfishman.github.io/ITensorChemistry.jl/dev)
[![Build Status](https://github.com/mtfishman/ITensorChemistry.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mtfishman/ITensorChemistry.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mtfishman/ITensorChemistry.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mtfishman/ITensorChemistry.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Installation

```julia
julia> using Pkg

julia> Pkg.add(; url="https://github.com/mtfishman/ITensorChemistry.jl")
```

## Example usage

Run DMRG on a specified molecule in the molecular orbital basis:
```julia
using ITensors
using ITensorChemistry

molecule = "H₂"
basis = "sto-3g"

@show molecule
@show basis

(; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian_and_state(; molecule, basis)

println("Basis set size = ", length(state))

s = siteinds("Electron", length(state); conserve_qns=true)
H = MPO(hamiltonian, s)
ψhf = MPS(s, state)

@show inner(ψhf, H, ψhf)
@show hartree_fock_energy

sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-6)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
dmrg_energy, ψ = dmrg(H, ψhf, sweeps)

@show dmrg_energy, hartree_fock_energy
@show dmrg_energy < hartree_fock_energy
```
