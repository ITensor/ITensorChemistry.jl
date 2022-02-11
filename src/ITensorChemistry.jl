module ITensorChemistry

using Fermi
using ITensors
using Suppressor

import Combinatorics: levicivita

include(joinpath("molecules", "molecule.jl"))
for file in readdir(joinpath(@__DIR__, "molecules", "molecules"); join=true)
  if endswith(file, ".jl")
    include(file)
  end
end
include("molecular_orbital_hamiltonian.jl")
include("jordan_wigner.jl")

export molecular_orbital_hamiltonian,pauli

end
