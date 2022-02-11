using ITensors
using ITensorChemistry
using Test

@testset "ITensorChemistry.jl" begin
  molecule = "H₂"
  basis = "sto-3g"
  (; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian(; molecule, basis)

  s = siteinds("Electron", length(state); conserve_qns=true)
  H = MPO(hamiltonian, s)
  ψhf = MPS(s, state)

  @test inner(ψhf, H, ψhf) ≈ hartree_fock_energy

  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-6)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  dmrg_energy, ψ = dmrg(H, ψhf, sweeps; outputlevel=0)

  @test dmrg_energy < hartree_fock_energy
end

@testset "electron-fermion mapping" begin
  molecule = "H₂"
  basis = "sto-3g"
  (; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian(; molecule, basis)
  se = siteinds("Electron", length(state); conserve_qns=false)
  He = MPO(hamiltonian, se)
  ψe = MPS(se, state)

  @test inner(ψe, He, ψe) ≈ hartree_fock_energy

  (; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian(; molecule, basis, sitetype = "Fermion")
  sf = siteinds("Fermion", length(state); conserve_qns=false)
  Hf = MPO(hamiltonian, sf)
  ψf = MPS(sf, state)

  @test inner(ψf, Hf, ψf) ≈ hartree_fock_energy
  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-6)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  Ee, ψe = dmrg(He, randomMPS(se), sweeps; outputlevel=0)
  
  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-6)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  Ef, ψf = dmrg(Hf, randomMPS(sf), sweeps; outputlevel=0)
  @test Ee ≈ Ef
end

@testset "fermion-qubit mapping" begin
  molecule = "H₂"
  basis = "sto-3g"
  (; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian(; molecule, basis, sitetype = "Fermion")
  
  sf = siteinds("Fermion", length(state); conserve_qns=false)
  Hf = MPO(hamiltonian, sf)
  ψf = MPS(sf, state)

  @test inner(ψf, Hf, ψf) ≈ hartree_fock_energy
  

  qubit_hamiltonian = jordanwigner(hamiltonian)
  sq = siteinds("Qubit", length(state))
  Hq = MPO(qubit_hamiltonian, sq)
  ψq = MPS(sq, state)

  @test inner(ψq, Hq, ψq) ≈ hartree_fock_energy
  
  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-6)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  Ef, _ = dmrg(Hf, randomMPS(sf), sweeps; outputlevel=0)
  
  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-6)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  Eq, _ = dmrg(Hq, randomMPS(sq), sweeps; outputlevel=0)
  @test Ef ≈ Eq
end

