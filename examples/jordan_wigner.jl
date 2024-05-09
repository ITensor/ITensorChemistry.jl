using ITensors, ITensorMPS
using ITensorChemistry

# generate Nitrogen molecule
molecule = Molecule("N₂")
basis = "sto-3g"

@show molecule

hf = molecular_orbital_hamiltonian(molecule; basis)
hamiltonian = hf.hamiltonian
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy

println("Number of orbitals = ", length(hartree_fock_state))
println("Number of fermionic operators = ", length(hamiltonian))
println("|HF⟩ = |", prod(string.(hartree_fock_state)), "⟩")

qubit_state = jordanwigner(hartree_fock_state)
qubit_hamiltonian = jordanwigner(hamiltonian)
println("Number of qubit operators = ", length(qubit_hamiltonian))
println("|HF⟩ = |", prod(string.(qubit_state)), "⟩")
