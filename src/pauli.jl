mutable struct Pauli
  c::Number
  σ::String
end

pauli(c::Number, σ::String) = Pauli(c, σ)
  pauli(σ::String) = Pauli(1.0, σ)

operator(p::Pauli)    = p.σ
operator(p::Pauli, n) = operator(p)[n]

coefficient(p::Pauli) = p.c
length(p::Pauli)      = length(p.σ)

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
  state > length(pauli) && return nothing
  σ = string(collect(operator(pauli))[state])
  return (σ, state+1)
end

function Base.resize!(pauli::Pauli, n::Int)
  pauli.σ = pauli.σ * "I"^(n-length(pauli))
  return pauli
end

# map paulis to integers (and back) using the following encoding 
# 'I' <=> 0  'X' <=> 1  'Y' <=> 2  'Z' <=> 3
integer_to_pauli(i::Int) = ['I','X','Y','Z'][i+1]

function pauli_to_integer(c::Char) 
  i = findfirst(x -> x == c, ['I','X','Y','Z'])
  isnothing(i) && error("operator not recognized")
  return i-1 
end

pauli_to_integer(pauli::Pauli) = 
  [pauli_to_integer(c) for c in collect(operator(pauli))]


function *(A₀::Pauli, α::Number)
  A = copy(A₀)
  A.c *= α
  return A
end

*(α::Number, p::Pauli) = *(p, α)

function pauliproduct(σa::Char, σb::Char)
  σa == σb && return (1.,'I')
  ('I' in σa) && return (1., σb)
  ('I' in σb) && return (1., σa)
  μa, μb = pauli_to_integer.((σa, σb))
  μc = first(filter(x -> x ∉ [μa,μb], [1,2,3]))
  return (im * levicivita([μa, μb, μc]), integer_to_pauli(μc))
end

# multiply two paulis
function *(A::Pauli, B::Pauli)
  c = coefficient(A) * coefficient(B)
  la = length(A)
  lb = length(B)
  L = max(la, lb)

  # add identities
  la < lb && resize!(A, L)
  la > lb && resize!(B, L)
  
  Π = [pauliproduct(operator(A, n), operator(B, n)) for n in 1:L]
  c′ = c * prod(first(x) for x in Π)
  σ′ = join(last(x) for x in Π)
  return pauli(c′, σ′)
end
