using ITensors, ITensorMPS
using ITensorChemistry
using ITensorParallel
using MPI

MPI.Init()

ITensors.BLAS.set_num_threads(1)
ITensors.Strided.disable_threads()
ITensorChemistry.Fermi.TBLIS.set_num_threads(1)

ITensors.enable_threaded_blocksparse()

@show Threads.nthreads()

molecule = Molecule("H₂O")
basis = "sto-3g"

# Number of parts to split the Hamiltonian terms into
nparts = 10

@show molecule
@show basis

println("\nRunning Hartree-Fock")
hf = @time molecular_orbital_hamiltonian(molecule; basis)
opsum = hf.hamiltonian
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy
println("Hartree-Fock complete")

nsites = length(hartree_fock_state)

println("Basis set size = ", nsites)

sites = siteinds("Electron", nsites; conserve_qns=true)

println("\nConstruct MPO")

nprocs = MPI.Comm_size(MPI.COMM_WORLD)
opsums = collect(Iterators.partition(partition(opsum, nparts), nparts ÷ nprocs))

process = MPI.Comm_rank(MPI.COMM_WORLD) + 1
nterms = length(opsums[process])
H = Vector{MPO}(undef, nterms)
Threads.@threads for j in 1:nterms
  H[j] = MPO(opsums[process][j], sites)
end
println("MPO constructed")

@show maxlinkdim.(H)

PH = MPISum(ThreadedProjMPOSum(H))

ψhf = MPS(sites, hartree_fock_state)

ITensors.inner(ψ, Hs::Vector, ϕ) = sum([inner(ψ, H, ϕ) for H in Hs])
@show inner(ψhf', H, ψhf)
@show hartree_fock_energy

dmrg_params = (nsweeps=10, maxdim=[100, 200], cutoff=1e-6, noise=[1e-6, 1e-7, 1e-8, 0.0])

println("\nRunning DMRG")
@show dmrg_params

e, ψ = dmrg(PH, ψhf; dmrg_params...)
println("DMRG complete")

MPI.Finalize()
