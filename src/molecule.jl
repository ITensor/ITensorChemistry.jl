struct MoleculeName{M} end

MoleculeName(mn::Symbol) = MoleculeName{mn}()
MoleculeName(mn::String) = MoleculeName(Symbol(mn))

molecule(mol::String) = molecule(MoleculeName(mol))
molecule(::MoleculeName{mn}) where {mn} = error("`molecule` not implemented for $mn")

macro MoleculeName_str(s)
    :(MoleculeName{$(Expr(:quote, Symbol(s)))})
end


struct Atom
  name::String
  coordinate::NTuple{3, Real}
end

# attributes
name(atom::Atom)       = atom.name
coordinate(atom::Atom) = atom.coordinate

# generate an atom with coordinates
atom(name::String, coordinate::NTuple) = 
  Atom(name, coordinate)

atom(name::String, x, y, z) = 
  Atom(name, (x, y, z))

# generate an atom at r = (0,0,0)
atom(name::String) = 
  Atom(name, (0.0, 0.0, 0.0))


struct Molecule
  atoms::Vector{Atom}
end

# empty molecule
molecule() = Molecule(Atom[])

# molecule composed by a set of atoms
molecule(atoms::Vector{Atom}) = 
  Molecule(atoms)

molecule(a⃗::Atom...) = 
  molecule([a⃗...])


# parse atomic data from vectors of tuple/vectors
molecule(atoms::Vector{<:Tuple})  = 
  molecule([atom(a[1], a[2:4]) for a in atoms])
molecule(atoms::Vector{<:Vector}) = 
  molecule(Tuple.(atoms))

Base.length(mol::Molecule) = 
  length(mol.atoms)
Base.getindex(mol::Molecule, i::Int) = 
  mol.atoms[i] 
Base.push!(mol::Molecule, atom::Atom) = 
  push!(mol.atoms, atom)

# parse data contained in `mol` to generate
# an input string encoding for HF in Fermi.jl
function parse_molecule(mol::Molecule)
  molstr = ""
  for a in 1:length(mol)
    atomname = name(mol[a])
    atomcoords = Float64.(coordinate(mol[a]))
    molstr *= atomname 
    for r in atomcoords
      molstr *= " " * string(r)
    end
    molstr *= "\n"
  end
  return molstr
end


