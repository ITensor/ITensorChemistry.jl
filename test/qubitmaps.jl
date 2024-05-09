using ITensors, ITensorMPS
using ITensorChemistry
using Test

@testset "Jordan Wigner mapping" begin
  molecule = Molecule("H₂")
  basis = "sto-3g"
  hf = molecular_orbital_hamiltonian(molecule; basis)
  hamiltonian = hf.hamiltonian
  hartree_fock_state = hf.hartree_fock_state
  hartree_fock_energy = hf.hartree_fock_energy

  electron_hilbert = siteinds("Electron", length(hartree_fock_state); conserve_qns=true)
  He = MPO(hamiltonian, electron_hilbert)
  ψ₀e = MPS(electron_hilbert, hartree_fock_state)
  @test inner(ψ₀e', He, ψ₀e) ≈ hartree_fock_energy
  electronE = copy(hartree_fock_energy)

  hf = molecular_orbital_hamiltonian(molecule; basis, sitetype="Fermion")
  hamiltonian = hf.hamiltonian
  hartree_fock_state = hf.hartree_fock_state
  hartree_fock_energy = hf.hartree_fock_energy

  fermion_hilbert = siteinds("Fermion", length(hartree_fock_state); conserve_qns=true)
  Hf = MPO(hamiltonian, fermion_hilbert)
  ψ₀f = MPS(fermion_hilbert, hartree_fock_state)
  @test inner(ψ₀f', Hf, ψ₀f) ≈ electronE

  hf = molecular_orbital_hamiltonian(molecule; basis, sitetype="Qubit")
  hamiltonian = hf.hamiltonian
  hartree_fock_state = hf.hartree_fock_state
  hartree_fock_energy = hf.hartree_fock_energy

  qubit_hilbert = siteinds("Qubit", length(hartree_fock_state); conserve_qns=true)
  Hq = MPO(hamiltonian, qubit_hilbert)
  ψ₀q = MPS(qubit_hilbert, hartree_fock_state)
  @test inner(ψ₀q', Hq, ψ₀q) ≈ electronE

  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-10)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

  Ee, _ = dmrg(He, ψ₀e, sweeps; outputlevel=0)
  Ef, _ = dmrg(Hf, ψ₀f, sweeps; outputlevel=0)
  Eq, _ = dmrg(Hq, ψ₀q, sweeps; outputlevel=0)
  @test Ee ≈ Ef
  @test Ee ≈ Eq
end
