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
  hf_energy = hf_wfn.energy
  Vnuc = aoints.molecule.Vnuc

  return (; hα, gα, nαocc, hf_energy, Vnuc)
end

function molecular_orbital_hamiltonian_and_state(; kwargs...)
  (; hα, gα, nαocc, hf_energy, Vnuc) = molecular_orbital_hamiltonian_coefficients(; kwargs...)
  nα = size(hα, 1)
  αocc = 1:nαocc

  # Representation of the second quantized quantum chemistry Hamiltonian.
  os = OpSum()
  os += Vnuc, "Id", 1
  for i in 1:nα, j in 1:nα
    if norm(hα[i, j]) > 1e-15
      os .+= hα[i, j], "Cdagup", i, "Cup", j
      os .+= hα[i, j], "Cdagdn", i, "Cdn", j
    end
  end
  for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα
    if norm(gα[i, j, k, l]) > 1e-15
      if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
        os .+= gα[i, j, k, l], "Cdagup", i, "Cdagup", j, "Cup", k, "Cup", l
        os .+= gα[i, j, k, l], "Cdagdn", i, "Cdagdn", j, "Cdn", k, "Cdn", l
      end
      os .+= gα[i, j, k, l], "Cdagup", i, "Cdagdn", j, "Cdn", k, "Cup", l
      os .+= gα[i, j, k, l], "Cdagdn", i, "Cdagup", j, "Cup", k, "Cdn", l
    end
  end

  αocc_state = [i in αocc ? 1 : 0 for i in 1:nα]
  occ_to_state = Dict([0 => 1, 1 => 4])
  init_state = [occ_to_state[n] for n in αocc_state]
  return (; hamiltonian=os, state=init_state, hartree_fock_energy=hf_energy)
end
