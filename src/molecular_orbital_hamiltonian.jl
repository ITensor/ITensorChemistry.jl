"""
Source: https://github.com/FermiQC/Fermi.jl/discussions/117

Credit to: Gustavo Aroeira (https://github.com/gustavojra)
"""
function molecular_orbital_hamiltonian_coefficients(;
  molecule,
  basis="sto-3g",
  diis=diis(Molecule(molecule)),
  oda=oda(Molecule(molecule))
)

  @suppress begin
    Fermi.Options.set("molstring", ITensorChemistry.molecule(Molecule(molecule)))
    Fermi.Options.set("basis", basis)
    Fermi.Options.set("diis", diis)
    Fermi.Options.set("oda", oda)
  end
  hf_wfn = @suppress Fermi.HartreeFock.RHF()
  aoints = @suppress Fermi.Integrals.IntegralHelper(orbitals=hf_wfn.orbitals)

  hα = @suppress aoints["T"] + aoints["V"]
  gα = @suppress aoints["ERI"]

  # Convert `g` from Chemist convention to Physicist convention
  gα = 0.5 * permutedims(gα, (3, 2, 1, 4))

  nαocc = hf_wfn.molecule.Nα
  hartree_fock_energy = hf_wfn.energy
  nuclear_energy = aoints.molecule.Vnuc

  return (; hα, gα, nαocc, hartree_fock_energy, nuclear_energy)
end


function _molecular_orbital_hamiltonian(hα, gα, nuclear_energy)
  # Representation of the second quantized quantum chemistry Hamiltonian.
  hamiltonian = OpSum()
  hamiltonian += nuclear_energy, "Id", 1
  nα = size(hα, 1)
  for i in 1:nα, j in 1:nα
    if norm(hα[i, j]) > 1e-15
      hamiltonian .+= hα[i, j], "c†↑", i, "c↑", j
      hamiltonian .+= hα[i, j], "c†↓", i, "c↓", j
    end
  end
  for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα
    if norm(gα[i, j, k, l]) > 1e-15
      if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
        hamiltonian .+= gα[i, j, k, l], "c†↑", i, "c†↑", j, "c↑", k, "c↑", l
        hamiltonian .+= gα[i, j, k, l], "c†↓", i, "c†↓", j, "c↓", k, "c↓", l
      end
      hamiltonian .+= gα[i, j, k, l], "c†↑", i, "c†↓", j, "c↓", k, "c↑", l
      hamiltonian .+= gα[i, j, k, l], "c†↓", i, "c†↑", j, "c↑", k, "c↓", l
    end
  end
  return hamiltonian
end

function proc(p, ptot, i)
  return p == mod1(i, ptot)
end

# i ≥ j
f(i, j) = (i - 1) * i ÷ 2 + j
function proc(p, ptot, i, j)
  if i ≤ j
    return p == mod1(f(j, i), ptot)
  end
  return p == mod1(f(i, j), ptot)
end

function proc_sorted(p, ptot, i, j, k, l)
  if j == k
    return proc(p, ptot, j)
  end
  return proc(p, ptot, i, j)
end

function proc(p, ptot, i, j, k, l)
  return proc_sorted(p, ptot, sort((i, j, k, l))...)
end

# From: https://arxiv.org/abs/2103.09976
function _molecular_orbital_hamiltonian(hα, gα, nuclear_energy, nsub_hamiltonians)
  # Representation of the second quantized quantum chemistry Hamiltonian.
  # Split over `nsub_hamiltonians` to be parallelized over.
  hamiltonian = [OpSum() for _ in 1:nsub_hamiltonians]
  hamiltonian[1] += nuclear_energy, "Id", 1
  nα = size(hα, 1)
  for p in 1:nsub_hamiltonians
    for i in 1:nα, j in 1:nα
      pij = 1 / 2 * (proc(p, nsub_hamiltonians, i) + proc(p, nsub_hamiltonians, j))
      if pij ≠ 0 && norm(hα[i, j]) > 1e-10
        hamiltonian[p] .+= pij * hα[i, j], "c†↑", i, "c↑", j
        hamiltonian[p] .+= pij * hα[i, j], "c†↓", i, "c↓", j
      end
    end
    for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα
      pijkl = proc(p, nsub_hamiltonians, i, j, k, l)
      if pijkl ≠ 0 && norm(gα[i, j, k, l]) > 1e-10
        if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
          hamiltonian[p] .+= gα[i, j, k, l], "c†↑", i, "c†↑", j, "c↑", k, "c↑", l
          hamiltonian[p] .+= gα[i, j, k, l], "c†↓", i, "c†↓", j, "c↓", k, "c↓", l
        end
        hamiltonian[p] .+= gα[i, j, k, l], "c†↑", i, "c†↓", j, "c↓", k, "c↑", l
        hamiltonian[p] .+= gα[i, j, k, l], "c†↓", i, "c†↑", j, "c↑", k, "c↓", l
      end
    end
  end
  return hamiltonian
end

function molecular_orbital_hamiltonian(nsub_hamiltonians=nothing; sitetype::String = "Electron", kwargs...)
  (; hα, gα, nαocc, hartree_fock_energy, nuclear_energy) = molecular_orbital_hamiltonian_coefficients(; kwargs...)

  if isnothing(nsub_hamiltonians)
    hamiltonian = _molecular_orbital_hamiltonian(hα, gα, nuclear_energy)
  else
    hamiltonian = _molecular_orbital_hamiltonian(hα, gα, nuclear_energy, nsub_hamiltonians)
  end

  nα = size(hα, 1)
  αocc_state = [i in 1:nαocc ? 1 : 0 for i in 1:nα]
  occ_to_state = Dict([0 => 1, 1 => 4])
  state = [occ_to_state[n] for n in αocc_state]
  
  if sitetype == "Fermion"
    hamiltonian = electron_to_fermion(hamiltonian)
    state = electron_to_fermion(state)
    return (; hamiltonian, state, hartree_fock_energy)
  end
  return (; hamiltonian, state, hartree_fock_energy)
end


"""
    electron_to_fermion(hamiltonian::OpSum)

Map an OpSum from spinfull to spinless fermions.
"""
function electron_to_fermion(hamiltonian::OpSum)
  fermion_hamiltonian = OpSum()
  # loop over MPOTerms
  for k in 1:length(hamiltonian)
    h = hamiltonian[k]
    c = ITensors.coef(h)
    sites = first.(ITensors.sites.(ITensors.ops(h)))
    O = ITensors.name.(ITensors.ops(h))
    
    if O == ["Id"]
      fermion_hamiltonian += c, "Id", 1
    else
      ops_and_sites = []
      # loop over each single-site operator
      for (j,o) in enumerate(O)
        # the fermion is placed at twice the site number + 1 if ↓
        fermionsite = 2 * sites[j] + Int(o[end] == '↓') - 1
        ops_and_sites = vcat(ops_and_sites, (String(strip(o, o[end])), fermionsite)...) 
      end
      fermion_hamiltonian += (c, ops_and_sites...)
    end
  end
  return fermion_hamiltonian
end

"""
    electron_to_fermion(electron_state::Vector)

Map a product state from spinfull to spinless fermions.
"""
function electron_to_fermion(electron_state::Vector)
  stmap = [[1,1],[1,2],[2,1],[2,2]]
  fermion_state = zeros(Int,2*length(electron_state))
  for k in 1:length(electron_state)
    fermion_state[2*k-1:2*k] = stmap[electron_state[k]]
  end
  return fermion_state
end

