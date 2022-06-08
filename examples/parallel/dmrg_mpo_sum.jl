using ITensors
using ITensorChemistry
using ITensorParallel

ITensors.BLAS.set_num_threads(1)
ITensors.Strided.disable_threads()
ITensors.disable_threaded_blocksparse()

@show Threads.nthreads()

molecule = "H₂O"
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

hamiltonians = partition(hamiltonian, Threads.nthreads())
H = Vector{MPO}(undef, Threads.nthreads())
@time Threads.@threads for n in 1:Threads.nthreads()
  H[Threads.threadid()] = MPO(hamiltonians[n], s)
end
println("MPO constructed")

@show maxlinkdim.(H)

ψhf = MPS(s, hartree_fock_state)

ITensors.inner(ψ, Hs::Vector, ϕ) = sum([inner(ψ, H, ϕ) for H in Hs])
@show inner(ψhf', H, ψhf)
@show hartree_fock_energy

dmrg_params = (nsweeps=10, maxdim=[100, 200], cutoff=1e-6, noise=[1e-6, 1e-7, 1e-8, 0.0])

println("\nRunning DMRG")
@show dmrg_params

e, ψ = dmrg(ThreadedProjMPOSum(H), ψhf; dmrg_params...)
println("DMRG complete")
