using ITensors
using ITensorChemistry
using Test

@testset "fermion-qubit mapping" begin
  mol = molecule("H₂")
  basis = "sto-3g"
  (; hamiltonian, state, hartree_fock_energy) =
    molecular_orbital_hamiltonian(mol; basis, 
                                       diis = false, 
                                       oda = false, 
                                       sitetype = "Fermion")
  
  sf = siteinds("Fermion", length(state); conserve_qns=true)
  Hf = MPO(hamiltonian, sf)
  ψf = MPS(sf, state)

  @test inner(ψf, Hf, ψf) ≈ hartree_fock_energy
  

  qubit_hamiltonian = jordanwigner(hamiltonian)
  sq = siteinds("Qubit", length(state); conserve_qns=true)
  Hq = MPO(qubit_hamiltonian, sq)
  ψq = MPS(sq, state)

  @test inner(ψq, Hq, ψq) ≈ hartree_fock_energy
  
  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-10)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  Ef, _ = dmrg(Hf, ψf, sweeps; outputlevel=0)
  
  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-10)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  Eq, _ = dmrg(Hq, ψq, sweeps; outputlevel=0)
  @test Ef ≈ Eq
end

