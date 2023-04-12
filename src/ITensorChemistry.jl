module ITensorChemistry

using ITensors
using PythonCall

import Combinatorics: levicivita

import Base: *, length, getindex, setindex!, iterate, copy, push!, resize!

include("molecule.jl")
include("molecules.jl")
include("molecular_orbital_hamiltonian.jl")
include("pauli.jl")
include("qubitmaps.jl")

const pyscf = PythonCall.pynew() # Init to NULL
function __init__()
  PythonCall.pycopy!(pyscf, pyimport("pyscf"))
  return nothing
end

export Atom, Molecule, molecular_orbital_hamiltonian, jordanwigner
end
