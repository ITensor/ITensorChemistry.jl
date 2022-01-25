using ITensors
using ITensorChemistry

ITensors.BLAS.set_num_threads(1)
ITensors.Strided.disable_threads()
ITensors.disable_threaded_blocksparse()

@show Threads.nthreads()

# This pirates ITensors.jl functions in DMRG
# to make them parallel over the sum of MPOs.
include("parallel_projmpo.jl")

molecule = "H₂O"
basis = "sto-3g"

@show molecule
@show basis

println("\nRunning Hartree-Fock")
(; hamiltonian, state, hartree_fock_energy) = @time molecular_orbital_hamiltonian(Threads.nthreads(); molecule, basis)
println("Hartree-Fock complete")

println("Basis set size = ", length(state))

s = siteinds("Electron", length(state); conserve_qns=true)

println("\nConstruct MPO")

H = Vector{MPO}(undef, Threads.nthreads())
@time Threads.@threads for n in 1:Threads.nthreads()
  H[Threads.threadid()] = MPO(hamiltonian[n], s)
end
println("MPO constructed")

@show maxlinkdim.(H)

ψhf = MPS(s, state)

ITensors.inner(ψ, Hs::Vector, ϕ) = sum([inner(ψ, H, ϕ) for H in Hs])
@show inner(ψhf, H, ψhf)
@show hartree_fock_energy

sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-6)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

println("\nRunning DMRG")
@show sweeps

e, ψ = dmrg(H, ψhf, sweeps)
println("DMRG complete")
