function jordanwigner(H::OpSum; threshold = 1e-10)
  Hq = OpSum()
  for k in 1:length(H)
    if ITensors.name.(ITensors.ops(H[k])) == ["Id"]
      Hq += ITensors.coef(H[k]), "Id", 1
    else
      Q = jordanwigner(H[k])
      for j in 1:length(Q)
        Hq += Q[j]
      end
    end
    ITensors.sortmergeterms!(Hq)
  end
  newH = OpSum()
  for k in 1:length(Hq)
    coeff = ITensors.coef(Hq[k])
    if norm(coeff) > threshold
      push!(newH, Hq[k])
    end
  end
  return newH
end

function jordanwigner(F::ITensors.MPOTerm)

  Fcoeff = ITensors.coef(F)
  Fops   = ITensors.ops(F)
  Fnames = ITensors.name.(Fops)
  Fsites = first.(ITensors.sites.(Fops))
  
  Fops == ["Id"] && return ITensors.MPOTerm(1.0, "Id", 1)
  
  oplist = []
  coefflist = []
  maxsite = maximum(Fsites)
  for k in 1:length(Fops)
    Fsite = Fsites[k]
    Fname = Fnames[k]
    fstring = pauli(repeat("Z", Fsite-1))
    X = pauli("X", Fsite)
    X = fstring * pauli(1/2, "X", Fsite)
    Y = fstring * pauli("Y", Fsite)
    Y = Fname[end] == '†' ? -im/2 * Y : im/2 * Y
    push!(oplist, (X,Y)) 
  end
  oplist = vec(collect(Base.Iterators.product(oplist...)))
  Qop = []
  energy_offset = 0.0
  for p in oplist
    O = reduce(*, p)
    coeff = (-1)^length(Fops) * Fcoeff * coefficient(O)
    op = operator(O)
    paulilocs = findall(x -> x ≠ 'I', collect(op))
    if !isempty(paulilocs)
      Qop = vcat(Qop, (coeff, vcat([vcat([string.(collect(op))[loc], loc]...) for loc in paulilocs]...)...)) 
    else
      energy_offset += coeff
    end
  end
  energy_offset ≠ 0.0 && return vcat(Qop, [(energy_offset, "Id", 1)]) 
  return Qop
end

