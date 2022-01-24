using ITensors
using ITensorChemistry
using Test

@testset "ITensorChemistry.jl" begin
  molecule = "H₂"
  basis = "sto-3g"
  (; hamiltonian, state, hartree_fock_energy) = molecular_orbital_hamiltonian_and_state(; molecule, basis)

  s = siteinds("Electron", length(state); conserve_qns=true)
  H = MPO(hamiltonian, s)
  ψhf = MPS(s, state)

  @test inner(ψhf, H, ψhf) ≈ hartree_fock_energy

  sweeps = Sweeps(10)
  setmaxdim!(sweeps, 100, 200)
  setcutoff!(sweeps, 1e-6)
  setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
  dmrg_energy, ψ = dmrg(H, ψhf, sweeps; output_level=0)

  @test dmrg_energy < hartree_fock_energy
end
