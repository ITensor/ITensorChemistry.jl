"""
Hs::Vector{MPO}
psi0::MPS
Hs = permute.(Hs, Ref((linkind, siteinds, linkind)))
PH = ExtendedProjMPOSum(Hs)
orthogonalize!(psi0, 1)
position!(PH, psi0, 1)
dmrg(PH, psi, sweeps; kwargs...)  # Keeping track of PHK allows one to restart without extra cost
"""
module ExtendedProjMPOs

using ITensors

import ITensors: product, position!, noiseterm, lproj, rproj

"""
Edge tensor with added on MPO tensor; srw 4/30/21
o--o--o-   
|  |  |  | 
o--o--o--o-
|  |  |  |
o--o--o- 
"""
mutable struct ExtendedProjMPO
  projmpo::ProjMPO
  EL::ITensor
  ER::ITensor
  ExtendedProjMPO(H::MPO) = new(ProjMPO(H), ITensor(), ITensor())
end

nsite(p::ExtendedProjMPO) = ITensors.nsite(p.projmpo)
Base.length(p::ExtendedProjMPO) = length(p.projmpo)
(p::ExtendedProjMPO)(v::ITensor) = product(p.projmpo, v)
Base.eltype(p::ExtendedProjMPO) = ITensors.eltype(p.projmpo)
Base.size(p::ExtendedProjMPO) = ITensors.size(p.projmpo)

lproj(p::ExtendedProjMPO) = ITensors.lproj(p.projmpo)
rproj(p::ExtendedProjMPO) = ITensors.rproj(p.projmpo)

function position!(p::ExtendedProjMPO, psi::MPS, pos::Int)
  ITensors.position!(p.projmpo, psi, pos)
  L, R = lproj(p), rproj(p)
  isnothing(L) && (L = 1.0)
  isnothing(R) && (R = 1.0)
  r = p.projmpo.rpos
  p.EL = L * p.projmpo.H[pos]
  p.ER = R * p.projmpo.H[r - 1]
  return p
end

function product(p::ExtendedProjMPO, v::ITensor)::ITensor
  Hv = p.EL * v * p.ER
  if order(Hv) != order(v)
    error(
      string(
        "The order of the ProjMPO-ITensor product P*v is not equal to the order of the ITensor v, ",
        "this is probably due to an index mismatch.\nCommon reasons for this error: \n",
        "(1) You are trying to multiply the ProjMPO with the $(nsite(p))-site wave-function at the wrong position.\n",
        "(2) `orthognalize!` was called, changing the MPS without updating the ProjMPO.\n\n",
        "P*v inds: $(inds(Hv)) \n\n",
        "v inds: $(inds(v))",
      ),
    )
  end
  return noprime(Hv)
end

function noiseterm(p::ExtendedProjMPO, phi::ITensor, ortho::String) # ITensors.noiseterm(p.projmpo,phi,ortho)
  if nsite(p) != 2
    error("noise term only defined for 2-site ProjMPO")
  end
  LR = (ortho == "left" ? p.EL : p.ER)
  U, S, V = svd(phi, commoninds(LR, phi)...; cutoff=1e-4)
  nt = LR * U * S
  return nt * dag(noprime(nt))
end

mutable struct ExtendedProjMPOSum
  pm::Vector{ExtendedProjMPO}
end

function ExtendedProjMPOSum(mpos::Vector{MPO})
  return ExtendedProjMPOSum([ExtendedProjMPO(M) for M in mpos])
end

ExtendedProjMPOSum(Ms::MPO...) = ExtendedProjMPOSum([Ms...])

nsite(P::ExtendedProjMPOSum) = nsite(P.pm[1])

Base.length(P::ExtendedProjMPOSum) = length(P.pm[1])

function product(P::ExtendedProjMPOSum, v::ITensor)::ITensor
  Pvs = fill(ITensor(), Threads.nthreads())
  Threads.@threads for Ppmn in P.pm
    Pvs[Threads.threadid()] += product(Ppmn, v)
  end
  return sum(Pvs)
end

function Base.eltype(P::ExtendedProjMPOSum)
  elT = eltype(P.pm[1])
  for n in 2:length(P.pm)
    elT = promote_type(elT, eltype(P.pm[n]))
  end
  return elT
end

(P::ExtendedProjMPOSum)(v::ITensor) = product(P, v)

Base.size(P::ExtendedProjMPOSum) = size(P.pm[1])

function position!(P::ExtendedProjMPOSum, psi::MPS, pos::Int)
  Threads.@threads for M in P.pm
    position!(M, psi, pos)
  end
  return P
end

function noiseterm(P::ExtendedProjMPOSum, phi::ITensor, dir::String)
  nts = fill(ITensor(), Threads.nthreads())
  Threads.@threads for pp in P.pm
    nts[Threads.threadid()] += noiseterm(pp, phi, dir)
  end
  return sum(nts)
end

export ExtendedProjMPO, ExtendedProjMPOSum

end
