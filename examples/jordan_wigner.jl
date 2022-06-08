using ITensors
using ITensorChemistry

# generate Nitrogen molecule
molecule = Molecule("N₂")
basis = "sto-3g"

@show molecule

(; hamiltonian, hartree_fock_state, hartree_fock_energy) = molecular_orbital_hamiltonian(;
  molecule, basis
)
println("Number of orbitals = ", length(hartree_fock_state))
println("Number of fermionic operators = ", length(hamiltonian))
println("|HF⟩ = |", prod(string.(hartree_fock_state)), "⟩")

qubit_state = jordanwigner(hartree_fock_state)
qubit_hamiltonian = jordanwigner(hamiltonian)
println("Number of qubit operators = ", length(qubit_hamiltonian))
println("|HF⟩ = |", prod(string.(qubit_state)), "⟩")
