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
      hamiltonian .+= hα[i, j], "Cdagup", i, "Cup", j
      hamiltonian .+= hα[i, j], "Cdagdn", i, "Cdn", j
    end
  end
  for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα
    if norm(gα[i, j, k, l]) > 1e-15
      if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
        hamiltonian .+= gα[i, j, k, l], "Cdagup", i, "Cdagup", j, "Cup", k, "Cup", l
        hamiltonian .+= gα[i, j, k, l], "Cdagdn", i, "Cdagdn", j, "Cdn", k, "Cdn", l
      end
      hamiltonian .+= gα[i, j, k, l], "Cdagup", i, "Cdagdn", j, "Cdn", k, "Cup", l
      hamiltonian .+= gα[i, j, k, l], "Cdagdn", i, "Cdagup", j, "Cup", k, "Cdn", l
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
        hamiltonian[p] .+= pij * hα[i, j], "Cdagup", i, "Cup", j
        hamiltonian[p] .+= pij * hα[i, j], "Cdagdn", i, "Cdn", j
      end
    end
    for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα
      pijkl = proc(p, nsub_hamiltonians, i, j, k, l)
      if pijkl ≠ 0 && norm(gα[i, j, k, l]) > 1e-10
        if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
          hamiltonian[p] .+= gα[i, j, k, l], "Cdagup", i, "Cdagup", j, "Cup", k, "Cup", l
          hamiltonian[p] .+= gα[i, j, k, l], "Cdagdn", i, "Cdagdn", j, "Cdn", k, "Cdn", l
        end
        hamiltonian[p] .+= gα[i, j, k, l], "Cdagup", i, "Cdagdn", j, "Cdn", k, "Cup", l
        hamiltonian[p] .+= gα[i, j, k, l], "Cdagdn", i, "Cdagup", j, "Cup", k, "Cdn", l
      end
    end
  end
  return hamiltonian
end

function molecular_orbital_hamiltonian(nsub_hamiltonians=nothing; kwargs...)
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
  return (; hamiltonian, state, hartree_fock_energy)
end
