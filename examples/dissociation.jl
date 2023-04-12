using ITensors
using ITensorChemistry
using Printf
using Plots

# bond distances
r⃗ = 0.3:0.03:3.0

# hilbert space
s = siteinds("Electron", 2; conserve_qns=true)

# dmrg params
sweeps = Sweeps(10)
setmaxdim!(sweeps, 10, 20, 30, 40, 50, 100)
setcutoff!(sweeps, 1e-8)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

function energy_at_bond(r)
  # define molecule geometry
  molecule = Molecule([("H", 0.0, 0.0, 0.0), ("H", r, 0.0, 0.0)])

  # build electronic hamiltonian and solve HF
  hf = molecular_orbital_hamiltonian(molecule; basis="sto-3g")
  hamiltonian = hf.hamiltonian
  hartree_fock_state = hf.hartree_fock_state
  hartree_fock_energy = hf.hartree_fock_energy

  H = MPO(hamiltonian, s)

  # initialize MPS to HF state
  ψhf = MPS(s, hartree_fock_state)

  # run dmrg
  dmrg_energy, _ = dmrg(H, ψhf, sweeps; outputlevel=0)
  return hartree_fock_energy, dmrg_energy
end

println("-"^40 * "\nr₁-r₂\t   Energy (HF)\t   Energy (DMRG)")
Ehf = [];
Edmrg = [];
for r in r⃗
  ehf, edmrg = energy_at_bond(r)
  @printf("%-11.3f%-16.8f%-16.8f\n", r, ehf, edmrg)
  push!(Ehf, ehf)
  push!(Edmrg, edmrg)
end

pargs = (
  size=(800, 400),
  dpi=1000,
  margin=5Plots.mm,
  title="Dissociation of Hydrogen Molecule",
  xlabel="bond distance (au)",
  ylabel="Ground State Energy",
  legend=:top,
  markersize=4,
  marker=:circle,
  linewidth=1,
  xguidefontsize=15,
  yguidefontsize=15,
  legendfontsize=12,
  titlefontsize=20,
  xtickfontsize=15,
  ytickfontsize=15,
)

p = plot(r⃗, Ehf; label="Hartree-Fock", pargs...)
p = plot(p, r⃗, Edmrg; label="DMRG", pargs...)
savefig(p, "dissociation.png")
p
