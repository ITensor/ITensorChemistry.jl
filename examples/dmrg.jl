using ITensors
using ITensorChemistry

molecule = "N₂"
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
e, ψ = dmrg(H, ψhf, sweeps)
