module ITensorChemistry

using Fermi
using ITensors
using Suppressor

import Combinatorics: levicivita

import Base: 
  sign, 
  *, 
  length, 
  getindex, 
  setindex!, 
  iterate, 
  copy, 
  push!, 
  resize!,
  product

include("molecule.jl")
include("molecules.jl")
include("molecular_orbital_hamiltonian.jl")
include("pauli.jl")
include("qubitmaps.jl")

export 
  Atom,
  Molecule,
  atom,
  molecule,
  molecular_orbital_hamiltonian,
  jordanwigner
end
