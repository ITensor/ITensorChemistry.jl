using ITensors, ITensorMPS
using ITensorChemistry
using Test

@testset "Molecules" begin
  mol = Molecule()
  @test length(mol) == 0

  a₁ = Atom("A", 0.0, 0.0, 0.0)
  a₂ = Atom("B", 0.0, 0.5, 1.0)
  a₃ = Atom("C", 0.0, -1.0, 2.0)

  mol = Molecule(a₁, a₂, a₃)
  @test length(mol) == 3

  @test mol[1] == a₁
  @test mol[2] == a₂
  @test mol[3] == a₃

  mol2 = Molecule([("A", 0.0, 0.0, 0.0), ("B", 0.0, 0.5, 1.0), ("C", 0.0, -1.0, 2.0)])
  for k in 1:length(mol)
    @test mol2[k] == mol[k]
  end

  @test length(Molecule("H₂")) == 2
  @test length(Molecule("N₂")) == 2
  @test length(Molecule("H₂O")) == 3
  @test length(Molecule("BeH₂")) == 3
  @test length(Molecule("NH₃")) == 4
  @test length(Molecule("CH₂O")) == 4
  @test length(Molecule("CH₄")) == 5
  @test length(Molecule("C₂H₅OH")) == 9
  @test length(Molecule("C₆H₆")) == 12
end
