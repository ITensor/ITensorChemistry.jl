using ITensors
using ITensorChemistry
using Test

@testset "molecules" begin
  mol = molecule()
  @test length(mol) == 0
  
  a₁ = atom("A", 0.0, 0.0,  0.0)
  a₂ = atom("B", 0.0, 0.5,  1.0)
  a₃ = atom("C", 0.0, -1.0, 2.0)
   
  mol = molecule(a₁, a₂, a₃)
  @test length(mol) == 3
  
  @test mol[1] == a₁
  @test mol[2] == a₂
  @test mol[3] == a₃
  
  
  mol2 = molecule([("A", 0.0, 0.0, 0.0),
                   ("B", 0.0, 0.5,  1.0),
                   ("C", 0.0, -1.0, 2.0)])
  for k in 1:length(mol)
    @test mol2[k] == mol[k]
  end

  @test length(molecule("H₂")) == 2
  @test length(molecule("N₂")) == 2
  @test length(molecule("H₂O")) == 3
  @test length(molecule("BeH₂")) == 3
  @test length(molecule("NH₃")) == 4
  @test length(molecule("CH₂O")) == 4
  @test length(molecule("CH₄")) == 5
  @test length(molecule("C₂H₅OH")) == 9
  @test length(molecule("C₆H₆")) == 12
end

