"""
    molecular_orbital_hamiltonian_coefficients
    molecule::String 
    basis="sto-3g"
    charge=0
    spin=0
"""
function molecular_orbital_hamiltonian_coefficients(
  molecule::String; basis="sto-3g", charge=0, spin=0
)
  # Input Check (TODO add UHF)
  if spin != 0
    error("ITensorChemistry can only handle spin=0 systems right now!")
  end

  mol = pyscf.gto.M(; atom=molecule, basis=basis, charge=charge, spin=spin, verbose=3)

  # Run HF
  mf = pyscf.scf.RHF(mol)
  mf.chkfile = ".tmpfile_pyscf_itensor" # Set this bc Python.tempfile can cause problems sometimes
  mf.kernel()
  println("RHF Energy (Ha): ", mf.e_tot)

  # Create shorthands for 1- and 2-body integrals in MO basis
  mo = pyconvert(Array, mf.mo_coeff)
  hcore_ao = pyconvert(Array, mf.get_hcore())

  n = size(mo, 1)
  one_body = mo' * hcore_ao * mo
  two_body = reshape(pyconvert(Array, mol.ao2mo(mf.mo_coeff; aosym=1)), n, n, n, n)

  # Collect data from HF calculation to return
  hα = one_body
  gα = 0.5 * permutedims(two_body, (3, 2, 1, 4))
  nαocc = pyconvert(Number, mol.nelec[1])
  nuclear_energy = pyconvert(Number, mf.energy_nuc())
  hartree_fock_energy = pyconvert(Number, mf.e_tot)

  # println(pyconvert(String, mf.chkfile))
  rm(pyconvert(String, mf.chkfile))

  return (; hα, gα, nαocc, hartree_fock_energy, nuclear_energy)
end

function molecular_orbital_hamiltonian_coefficients(molecule::Molecule; kwargs...)
  return molecular_orbital_hamiltonian_coefficients(
    xyz_string(Molecule(molecule)); kwargs...
  )
end

function _molecular_orbital_hamiltonian(hα, gα, nuclear_energy; atol=1e-15)
  # Representation of the second quantized quantum chemistry Hamiltonian.
  hamiltonian = OpSum()
  add!(hamiltonian, nuclear_energy, "Id", 1)
  nα = size(hα, 1)
  for i in 1:nα, j in 1:nα
    if norm(hα[i, j]) > atol
      add!(hamiltonian, hα[i, j], "c†↑", i, "c↑", j)
      add!(hamiltonian, hα[i, j], "c†↓", i, "c↓", j)
    end
  end
  for i in 1:nα, j in 1:nα, k in 1:nα, l in 1:nα
    if norm(gα[i, j, k, l]) > atol
      if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
        add!(hamiltonian, gα[i, j, k, l], "c†↑", i, "c†↑", j, "c↑", k, "c↑", l)
        add!(hamiltonian, gα[i, j, k, l], "c†↓", i, "c†↓", j, "c↓", k, "c↓", l)
      end
      add!(hamiltonian, gα[i, j, k, l], "c†↑", i, "c†↓", j, "c↓", k, "c↑", l)
      add!(hamiltonian, gα[i, j, k, l], "c†↓", i, "c†↑", j, "c↑", k, "c↓", l)
    end
  end
  return hamiltonian
end

function molecular_orbital_hamiltonian(molecule; sitetype::String="Electron", kwargs...)
  res = molecular_orbital_hamiltonian_coefficients(molecule; kwargs...)
  hα, gα, nαocc, hartree_fock_energy, nuclear_energy = res.hα,
  res.gα, res.nαocc, res.hartree_fock_energy,
  res.nuclear_energy

  hamiltonian = _molecular_orbital_hamiltonian(hα, gα, nuclear_energy)

  nα = size(hα, 1)
  αocc_state = [i in 1:nαocc ? 1 : 0 for i in 1:nα]
  occ_to_state = Dict([0 => 1, 1 => 4])
  hartree_fock_state = [occ_to_state[n] for n in αocc_state]

  if sitetype == "Fermion"
    hartree_fock_state = electron_to_fermion(hartree_fock_state)
    hamiltonian = electron_to_fermion(hamiltonian)
  elseif sitetype == "Qubit"
    hartree_fock_state = jordanwigner(hartree_fock_state)
    hamiltonian = jordanwigner(hamiltonian)
  end
  return (; hamiltonian, hartree_fock_state, hartree_fock_energy)
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
    c = ITensors.coefficient(h)
    sites = first.(ITensors.sites.(ITensors.terms(h)))
    O = ITensors.name.(ITensors.terms(h))

    if O == ["Id"]
      fermion_hamiltonian += c, "Id", 1
    else
      ops_and_sites = []
      # loop over each single-site operator
      for (j, o) in enumerate(O)
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
  stmap = [[1, 1], [1, 2], [2, 1], [2, 2]]
  fermion_state = zeros(Int, 2 * length(electron_state))
  for k in 1:length(electron_state)
    fermion_state[(2 * k - 1):(2 * k)] = stmap[electron_state[k]]
  end
  return fermion_state
end
