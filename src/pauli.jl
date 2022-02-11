mutable struct Pauli
  c::Number
  σ::String
end

pauli(c::Number, σ::String) = Pauli(c, σ)
  pauli(σ::String) = Pauli(1.0, σ)

function pauli(c::Number, σ::String, site::Int) 
  length(σ) == 1 && return Pauli(c, repeat("I", site-1) * σ)
  error("")
end

pauli(σ::String, site::Int) = 
  pauli(1.0, σ, site)


operator(p::Pauli) = p.σ
coefficient(p::Pauli) = p.c
length(p::Pauli) = length(p.σ)

copy(pauli::Pauli) =
  Pauli(coefficient(pauli), operator(pauli))


getindex(p::Pauli, i::Int) = string(p.σ[i])
setindex!(pauli::Pauli, s::String, i::Int) = 
  setindex!(pauli, only(s), i)

function setindex!(pauli::Pauli, s::Char, i::Int)
  σ = collect(operator(pauli))
  σ[i] = s
  pauli.σ = join(σ)
  return pauli
end


function Base.iterate(pauli::Pauli, state=1)
  if state > length(pauli)
    return nothing
  else
    σ = string(collect(operator(pauli))[state])
    return (σ, state+1)
  end
end

"""
Pad with identities to match another pauli
"""
function enlarge!(pauli::Pauli, N::Int)
  pauli.σ = pauli.σ * repeat("I",N-length(pauli))  
  return pauli
end

# some encodings
integer_to_pauli(i::Int) = 
  ['I','X','Y','Z'][i+1]

function pauli_to_integer(c::Char) 
  i = findfirst(x -> x == c, ['I','X','Y','Z'])
  isnothing(i) && error("pauli not recognized")
  return i-1
end

pauli_to_integer(pauli::Pauli) = 
  [pauli_to_integer(c) for c in collect(operator(pauli))]

function *(A₀::Pauli, α::Number)
  A = copy(A₀)
  A.c *= α
  return A
end

*(α::Number, p::Pauli) = *(p,α)
# multiply two paulis
function *(A::Pauli, B::Pauli)
  c = coefficient(A) * coefficient(B)
  la = length(A)
  lb = length(B)
 
  # add identities
  la < lb && enlarge!(A, lb)
  la > lb && enlarge!(B, la)
  
  ia = pauli_to_integer(A)
  ib = pauli_to_integer(B)
  
  # new pauli
  σ = []
  for (j,μa) in enumerate(ia)
    μb = ib[j]
    # if both or either are 0 (Id), return as is, else
    # return the missing pauli from {X,Y,Z}
    μa == μb 
    μc = μa == μb ? 0  : 
         (0 ∈ μa) ? μb : 
         (0 ∈ μb) ? μa :
         first(filter(x -> x ∉ [μa,μb], [1,2,3]))
    # get the coefficient from the pauli permutation
    c = (μa == μb || μc == μa || μc == μb) ? c : c * im * levicivita([μa, μb, μc])
    push!(σ, integer_to_pauli(μc))
  end
  return pauli(c, join(σ))
end

