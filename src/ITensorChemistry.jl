module ITensorChemistry

using Fermi
using ITensors
using Suppressor

import Combinatorics: levicivita
import Base: sign, *, length, getindex, setindex!, iterate, copy

include(joinpath("molecules", "molecule.jl"))
for file in readdir(joinpath(@__DIR__, "molecules", "molecules"); join=true)
  if endswith(file, ".jl")
    include(file)
  end
end
include("molecular_orbital_hamiltonian.jl")
include("pauli.jl")
include("jordan_wigner.jl")

export 
  molecular_orbital_hamiltonian,
  pauli,
  jordanwigner

end
