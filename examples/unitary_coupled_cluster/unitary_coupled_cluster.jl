using ITensors
using ITensorChemistry
#using ChainRulesCore
#using OptimKit
#using Zygote

molecule = "H₂O"
basis = "sto-3g"

@show molecule
@show basis

# Split the Hamiltonian into `nparts` sub-Hamiltonians
nparts = 2

println("\nRunning Hartree-Fock")
(; hamiltonian, state, hartree_fock_energy) = @time molecular_orbital_hamiltonian(nparts; molecule, basis)
println("Hartree-Fock complete")

println("Basis set size = ", length(state))

s = siteinds("Electron", length(state); conserve_qns=true)

println("\nConstruct MPO")

H = Vector{MPO}(undef, nparts)
@time for n in 1:nparts
  H[n] = MPO(hamiltonian[n], s)
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

# Unitary coupled cluster (UCC)
T = []
nα = length(ψ)
for i in 1:nα, j in (i + 1):nα
  Tn = OpSum()
  Tn += "c†↑", i, "c↑", j
  Tn .-= "c†↑", j, "c↑", i
  push!(T, Tn)

  Tn = OpSum()
  Tn += "c†↓", i, "c↓", j
  Tn -= "c†↓", j, "c↓", i
  push!(T, Tn)
end
for i in 1:nα, j in (i + 1):nα, k in (j + 1):nα, l in (k + 1):nα
  Tn = OpSum()
  Tn += "c†↑", i, "c†↑", j, "c↑", k, "c↑", l
  Tn -= "c†↑", l, "c†↑", k, "c↑", j, "c↑", i
  push!(T, Tn)

  Tn = OpSum()
  Tn += "c†↓", i, "c†↓", j, "c↓", k, "c↓", l
  Tn -= "c†↓", l, "c†↓", k, "c↓", j, "c↓", i
  push!(T, Tn)

  Tn = OpSum()
  Tn += "c†↑", i, "c†↓", j, "c↓", k, "c↑", l
  Tn -= "c†↑", l, "c†↓", k, "c↓", j, "c↑", i
  push!(T, Tn)

  Tn = OpSum()
  Tn += "c†↓", i, "c†↑", j, "c↑", k, "c↓", l
  Tn -= "c†↓", l, "c†↑", k, "c↑", j, "c↓", i
  push!(T, Tn)
end

nt = length(T)
t = zeros(nt)
U = [exp(t[n] * T[n]) for n in 1:nt]
#∂ₜU = [T[n] * exp(t[n] * T[n]) for n in 1:nt]

