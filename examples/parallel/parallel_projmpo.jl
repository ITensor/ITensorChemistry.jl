import ITensors: product, position!, noiseterm

function product(P::ProjMPOSum, v::ITensor)::ITensor
  Pvs = fill(ITensor(), Threads.nthreads())
  Threads.@threads for M in P.pm
    Pvs[Threads.threadid()] += product(M, v)
  end
  return sum(Pvs)
end

function position!(P::ProjMPOSum, psi::MPS, pos::Int)
  Threads.@threads for M in P.pm
    position!(M, psi, pos)
  end
  return P
end

function noiseterm(P::ProjMPOSum, phi::ITensor, dir::String)
  nts = fill(ITensor(), Threads.nthreads())
  Threads.@threads for M in P.pm
    nts[Threads.threadid()] += noiseterm(M, phi, dir)
  end
  return sum(nts)
end
