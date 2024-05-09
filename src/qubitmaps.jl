using ITensorMPS: ITensorMPS
"""
    jordanwigner(H::OpSum; threshold = 1e-15)

JORDAN-WIGNER TRANSFORMATION

aⱼ  = - (∏ᵢ₌₁ʲ  Zᵢ) (Xᵢ + iYᵢ)/2
aⱼ† = - (∏ᵢ₌₁ʲ  Zᵢ) (Xᵢ - iYᵢ)/2
"""
function jordanwigner(H::OpSum; atol=1e-20)
  # check if H needs to be transformed into Fermions
  Hf = _is_electronic(H) ? electron_to_fermion(H) : H
  Hq = OpSum()
  for k in 1:length(Hf)
    if ITensors.name.(ITensors.terms(Hf[k])) == ["Id"]
      Hq += ITensors.coefficient(Hf[k]), "Id", 1
    else
      jwops = jordanwigner(Hf[k])
      for j in 1:length(jwops)
        Hq += jwops[j]
      end
    end
  end
  Hq = ITensorMPS.sortmergeterms(Hq)
  prunedHq = OpSum()
  for k in 1:length(Hq)
    if norm(ITensors.coefficient(Hq[k])) > atol
      prunedHq += Hq[k]
    end
  end
  return prunedHq
end

function jordanwigner(F::Scaled{C,Prod{Op}}) where {C}
  Fcoeff = ITensors.coefficient(F)
  Fops = ITensors.terms(F)
  Fnames = ITensors.name.(Fops)
  Fsites = first.(ITensors.sites.(Fops))

  # return identity MPOTerm
  Fops == ["Id"] && return 1.0 * Prod{Op}() * ("I", 1)

  oplist = []
  coefflist = []
  for k in 1:length(Fops)
    Fsite = Fsites[k]
    Fname = Fnames[k]

    β = im / 2 * (-1)^Int(Fname[end] == '†')
    X = pauli(1 / 2, "Z"^(Fsite - 1) * "X")
    Y = pauli(β, "Z"^(Fsite - 1) * "Y")
    oplist = vcat(oplist, (X, Y))
  end
  Qop = []
  energy_offset = 0.0
  for p in Iterators.product(oplist...)
    P = reduce(*, p)
    coeff = (-1)^length(Fops) * Fcoeff * coefficient(P)
    O = string.(operator(P))
    paulilocs = findall(x -> x ≠ 'I', collect(O))
    if !isempty(paulilocs)
      Qop = vcat(
        Qop, (coeff, Tuple(vcat([[string(O[loc]), loc] for loc in paulilocs]...))...)
      )
    else
      energy_offset += coeff
    end
  end
  energy_offset ≠ 0.0 && return vcat(Qop, [(energy_offset, "Id", 1)])

  return Qop
end

jordanwigner(electron_state::Vector) = electron_to_fermion(electron_state)

function _is_electronic(H::OpSum)
  oplist = join(vcat([ITensors.name.(ITensors.terms(H[k])) for k in 1:length(H)]...))
  return ('↑' in oplist || '↓' in oplist)
end
