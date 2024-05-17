using PyCall
using ITensors, ITensorMPS

pyscf = pyimport("pyscf")
fci = pyimport("pyscf.fci")

#
# Create Simple Molecule
#
# mol = pyscf.gto.M(atom = "N 0 0 0; N 0 0 1;", basis = "sto3g", verbose = 3)
# mol = pyscf.gto.M(atom = "o 0 0 0; o 0 0 1;", basis = "sto3g", verbose = 3)
mol = pyscf.gto.M(; atom="""C      0.00000    0.00000    0.00000
  H      0.00000    0.00000    1.08900
  H      1.02672    0.00000   -0.36300
  H     -0.51336   -0.88916   -0.36300
  H     -0.51336    0.88916   -0.36300""", basis="sto3g", verbose=3)

# Run HF
mf = pyscf.scf.RHF(mol).run()
println("RHF Energy (Ha): ", mf.e_tot)

# Create shorthands for 1- and 2-body integrals in MO basis
mo = mf.mo_coeff
n = size(mo, 1)
one_body = mo' * mf.get_hcore() * mo
two_body = reshape(mol.ao2mo(mf.mo_coeff; aosym=1), n, n, n, n)

# FCI (i.e. exact diagonalization)
cisolver = fci.FCI(mf)
cisolver.kernel()
println("FCI Energy (Ha): ", cisolver.e_tot)

#
# Setup for MPS Calculation
#
t = one_body
V = 0.5 * permutedims(two_body, (3, 2, 1, 4))
n_occ = mf.mo_occ
e_nuclear = mf.energy_nuc()

# Create operators and MPO
function chemistry_hamiltonian(; t, V, e_nuclear)
  hamiltonian = OpSum()
  hamiltonian += e_nuclear, "Id", 1
  for i in 1:n, j in 1:n
    hamiltonian += t[i, j], "Cdagup", i, "Cup", j
    hamiltonian += t[i, j], "Cdagdn", i, "Cdn", j
  end
  for i in 1:n, j in 1:n, k in 1:n, l in 1:n
    if norm(V[i, j, k, l]) > 1e-15
      if (i ≠ j) && (k ≠ l) # Otherwise the terms are exactly zero
        hamiltonian += V[i, j, k, l], "Cdagup", i, "Cdagup", j, "Cup", k, "Cup", l
        hamiltonian += V[i, j, k, l], "Cdagdn", i, "Cdagdn", j, "Cdn", k, "Cdn", l
      end
      hamiltonian += V[i, j, k, l], "Cdagup", i, "Cdagdn", j, "Cdn", k, "Cup", l
      hamiltonian += V[i, j, k, l], "Cdagdn", i, "Cdagup", j, "Cup", k, "Cdn", l
    end
  end
  return hamiltonian
end

hamiltonian = chemistry_hamiltonian(; t, V, e_nuclear)
s = siteinds("Electron", n; conserve_qns=true)
H = MPO(hamiltonian, s)

#
# DMRG
#
#
# Check that we can make an MPS properly
occupation_to_state(n) = n == 0 ? 1 : (n == 2 ? 4 : error("Occupation is $n"))
ψmf = MPS(s, occupation_to_state.(n_occ))
e_mf_mps = inner(ψmf', H, ψmf)
println("Energy Error from MF MPS (Ha) ", abs(e_mf_mps - mf.e_tot))

# Initialize our MPS
ψ0 = random_mps(s, occupation_to_state.(n_occ); linkdims=40)
@show inner(ψ0', H, ψ0)

# Run DMRG
dmrg_params = (nsweeps=10, maxdim=[100, 300], cutoff=1e-7, noise=[1e-6, 1e-7, 1e-8, 0.0])
e, ψ = dmrg(H, ψ0; dmrg_params...)
println("DMRG Error ", abs(e - cisolver.e_tot))
