using Test

@testset "ITensorChemistry.jl" begin
  @testset "$filename" for filename in (
    "molecules.jl", "pauli.jl", "molecular_orbital_hamiltonian.jl", "qubitmaps.jl"
  )
    println("Running $filename")
    include(filename)
  end
end
