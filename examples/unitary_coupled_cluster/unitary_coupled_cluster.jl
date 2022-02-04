using ITensors
using ITensorChemistry
using OptimKit
using Random
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

function ITensors.op(::OpName"ucc_1e", st::SiteType"Electron", s1::Index, s2::Index; t)
  T = zeros(4^2, 4^2)
  for σ in ["↑", "↓"]
    # Note: kron has opposite conventions from ITensor
    #T += c†σ cσ
    T += kron(op("c$σ", st), op("c†$σ", st))
    #T -= cσ c†σ
    T -= kron(op("c†$σ", st), op("c$σ", st))
  end
  return itensor(exp(t * T), s1', s2', dag(s1), dag(s2))
end

function ITensors.op(::OpName"ucc_2e", st::SiteType"Electron", s1::Index, s2::Index, s3::Index, s4::Index; t)
  T = zeros(4^4, 4^4)
  for σ1 in ["↑", "↓"], σ2 in ["↑", "↓"]
    # Note: kron has opposite conventions from ITensor
    #T += c†σ c†σ cσ cσ
    T += kron(op("c$σ1", st), op("c$σ2", st), op("c†$σ2", st), op("c†$σ1", st))
    #T -= cσ cσ c†σ c†σ
    T -= kron(op("c†$σ1", st), op("c†$σ2", st), op("c$σ2", st), op("c$σ1", st))
  end
  return itensor(exp(t * T), s1', s2', s3', s4', dag(s1), dag(s2), dag(s3), dag(s4))
end

# Get the number of parameters
nα = length(s)
sites_1e = [(i, j) for i in 1:nα, j in 1:nα if i < j]
sites_2e = [(i, j, k, l) for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα if i < j < k < l]

function ucc_circuit(s, t)
  nt_1e = length(sites_1e)
  nt_2e = length(sites_2e)
  t_1e = t[1:nt_1e]
  t_2e = t[(nt_1e + 1):end]

  U_1e = [(n = sites_1e[i]; op("ucc_1e", s[n[1]], s[n[2]]; t=t_1e[i])) for i in eachindex(sites_1e)]
  U_2e = [(n = sites_2e[i]; op("ucc_2e", s[n[1]], s[n[2]], s[n[3]], s[n[4]]; t=t_2e[i])) for i in eachindex(sites_2e)]
  return [U_1e; U_2e]
end

function loss(t)
  U = ucc_circuit(s, t)
  Uψhf = apply(U, ψhf; cutoff=1e-6)
  return inner(Uψhf, H, Uψhf)
end

nt_1e = length(sites_1e)
nt_2e = length(sites_2e)

# Start at HF GS
# When trying to optimize starting from this state,
# for some reason the paramaters go to NaN.
# Need to debug.
#t0_1e = zeros(nt_1e)
#t0_2e = zeros(nt_2e)

Random.seed!(1234)
t0_1e = randn(nt_1e)
t0_2e = randn(nt_2e)

t0 = [t0_1e; t0_2e]

U0 = ucc_circuit(s, t0)

U0ψhf = apply(U0, ψhf)
@show inner(U0ψhf, H, U0ψhf)

@show loss(t0)

loss_∇loss(t) = (loss(t), convert(Vector, loss'(t)))
algorithm = LBFGS(; gradtol=1e-3, verbosity=2)
tₒₚₜ, lossₒₚₜ, ∇lossₒₚₜ, numfg, normgradhistory = optimize(loss_∇loss, t0, algorithm)

@show loss(tₒₚₜ)
