using ITensors, ITensorMPS
using ITensorChemistry
using Test

@testset "pauli" begin
  P = ITensorChemistry.pauli("X")
  @test length(P) == 1
  @test ITensorChemistry.coefficient(P) ≈ 1.0
  @test ITensorChemistry.operator(P) == 'X'

  P = ITensorChemistry.pauli(im, "X")
  @test ITensorChemistry.coefficient(P) ≈ im

  P = ITensorChemistry.pauli("XYZ")
  @test length(P) == 3
  @test ITensorChemistry.operator(P) == "XYZ"
  @test P[1] == 'X'
  @test P[2] == 'Y'
  @test P[3] == 'Z'

  q = ['X', 'Y', 'Z']
  for (i, p) in enumerate(P)
    @test p == q[i]
  end

  P = ITensorChemistry.pauli("XYZ")
  ITensorChemistry.resize!(P, 10)
  for i in 4:10
    @test P[i] == 'I'
  end
end

@testset "pauli product" begin
  X = ITensorChemistry.pauli("X")
  Y = ITensorChemistry.pauli("Y")
  Z = ITensorChemistry.pauli("Z")

  P = X * X
  @test ITensorChemistry.coefficient(P) ≈ 1
  @test ITensorChemistry.operator(P) == 'I'
  P = Y * Y
  @test ITensorChemistry.coefficient(P) ≈ 1
  @test ITensorChemistry.operator(P) == 'I'
  P = Z * Z
  @test ITensorChemistry.coefficient(P) ≈ 1
  @test ITensorChemistry.operator(P) == 'I'

  P = X * Y
  @test ITensorChemistry.coefficient(P) ≈ im
  @test ITensorChemistry.operator(P) == 'Z'
  P = Y * X
  @test ITensorChemistry.coefficient(P) ≈ -im
  @test ITensorChemistry.operator(P) == 'Z'

  P = Y * Z
  @test ITensorChemistry.coefficient(P) ≈ im
  @test ITensorChemistry.operator(P) == 'X'
  P = Z * Y
  @test ITensorChemistry.coefficient(P) ≈ -im
  @test ITensorChemistry.operator(P) == 'X'

  P = Z * X
  @test ITensorChemistry.coefficient(P) ≈ im
  @test ITensorChemistry.operator(P) == 'Y'
  P = X * Z
  @test ITensorChemistry.coefficient(P) ≈ -im
  @test ITensorChemistry.operator(P) == 'Y'

  A = ITensorChemistry.pauli(0.5, "YZXX")
  B = ITensorChemistry.pauli("XIZX")
  C = A * B

  @test ITensorChemistry.operator(C) == "ZZYI"
  @test ITensorChemistry.coefficient(C) ≈ -0.5
end
