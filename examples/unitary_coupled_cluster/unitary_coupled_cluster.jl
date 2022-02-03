using ITensors
using ITensorChemistry
using OptimKit
using Zygote

ITensors.inner(ψ, Hs::Vector, ϕ) = sum([inner(ψ, H, ϕ) for H in Hs])

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

@show inner(ψhf, H, ψhf)
@show hartree_fock_energy

sweeps = Sweeps(4)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-6)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

println("\nRunning DMRG")
@show sweeps

e, ψ = dmrg(H, ψhf, sweeps)
println("DMRG complete")

function ITensors.op(::OpName"ucc_1e", ::SiteType"Electron", s1::Index, s2::Index)
  T = ITensor()
  for σ in ["↑", "↓"]
    T += op("c†$σ", s1) * op("c$σ", s2)
    T -= op("c†$σ", s2) * op("c$σ", s1)
  end
  return T
end

function ITensors.op(::OpName"ucc_2e", ::SiteType"Electron", s1::Index, s2::Index, s3::Index, s4::Index)
  T = ITensor()
  for σ1 in ["↑", "↓"], σ2 in ["↑", "↓"]
    T += op("c†$σ1", s1) * op("c†$σ2", s2) * op("c$σ2", s3) * op("c$σ1", s4)
    T -= op("c†$σ1", s4) * op("c†$σ2", s3) * op("c$σ2", s2) * op("c$σ1", s1)
  end
  return T
end

nα = length(ψ)

T_1e = ITensor[op("ucc_1e", s[i], s[j]) for i in 1:nα, j in 1:nα if i < j]
T_2e = ITensor[op("ucc_2e", s[i], s[j], s[k], s[l]) for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα if i < j < k < l]
T = [T_1e; T_2e]

nt = length(T)
t0 = zeros(nt)
U = [exp(t0[n] * T[n]) for n in 1:nt]

@show length(U)

Uψhf = apply(U, ψhf)
@show inner(Uψhf, H, Uψhf)

function ucc_circuit(nα, t)
  # This is a more compact alternative, but Zygote fails to differentiate it
  #T_1e = ITensor[op("ucc_1e", s[i], s[j]) for i in 1:nα, j in 1:nα if i < j]
  T_1e = ITensor[]
  for i in 1:nα, j in (i + 1):nα
    T_1e = [T_1e; op("ucc_1e", s[i], s[j])]
  end
  T_2e = ITensor[]
  for i in 1:nα, j in (i + 1):nα, k in (j + 1):nα, l in (k + 1):nα
    T_2e = [T_2e; op("ucc_2e", s[i], s[j], s[k], s[l])]
  end
  T = [T_1e; T_2e]
  return [exp(t[n] * T[n]) for n in 1:length(t)]
end

function loss(t)
  U = ucc_circuit(nα, t)
  Uψhf = apply(U, ψhf)
  return inner(Uψhf, H, Uψhf)
end

@show loss(t0)

loss_∇loss(t) = (loss(t), convert(Vector, loss'(t)))
algorithm = LBFGS(; gradtol=1e-3, verbosity=2)
tₒₚₜ, lossₒₚₜ, ∇lossₒₚₜ, numfg, normgradhistory = optimize(loss_∇loss, t0, algorithm)

@show loss(tₒₚₜ)
