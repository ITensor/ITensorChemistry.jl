using ITensors, ITensorMPS
using ITensorChemistry

molecule = Molecule("H₂O")
basis = "sto-3g"

@show molecule
@show basis

println("\nRunning Hartree-Fock")
hf = @time molecular_orbital_hamiltonian(molecule; basis)
hamiltonian = hf.hamiltonian
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy

println("Hartree-Fock complete")

println("Basis set size = ", length(hartree_fock_state))

s = siteinds("Electron", length(hartree_fock_state); conserve_qns=true)

println("\nConstruct MPO")

H = @time MPO(hamiltonian, s)
println("MPO constructed")

@show maxlinkdim(H)

ψhf = MPS(s, hartree_fock_state)

@show inner(ψhf', H, ψhf)
@show hartree_fock_energy

dmrg_params = (nsweeps=10, maxdim=[100, 200], cutoff=1e-6, noise=[1e-6, 1e-7, 1e-8, 0.0])

println("\nRunning DMRG")
@show dmrg_params
e, ψ = dmrg(H, ψhf; dmrg_params...)
println("DMRG complete")
