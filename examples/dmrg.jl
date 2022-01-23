using ITensors
using ITensorChemistry

molecule = "N₂"

(; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian_and_state(; molecule)

s = siteinds("Electron", length(state); conserve_qns=true)
H = MPO(hamiltonian, s)
ψhf = MPS(s, state)

@show inner(ψhf, H, ψhf)
@show hartree_fock_energy

# TODO: add DMRG
