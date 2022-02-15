struct MoleculeName{M} end

MoleculeName(mn::Symbol) = MoleculeName{mn}()
MoleculeName(mn::String) = MoleculeName(Symbol(mn))

macro MoleculeName_str(s)
    :(MoleculeName{$(Expr(:quote, Symbol(s)))})
end


struct Atom
  name::String
  coordinates::NTuple{3, Real}
end

# attributes
name(atom::Atom)       = atom.name
coordinates(atom::Atom) = atom.coordinates

Atom(name::String, x, y, z) = 
  Atom(name, (x, y, z))

# generate an atom at r = (0,0,0)
Atom(name::String) = 
  Atom(name, (0.0, 0.0, 0.0))

struct Molecule
  atoms::Vector{Atom}
end

atoms(molecule::Molecule) = 
  molecule.atoms

# empty molecule
Molecule() = Molecule(Atom[])

Molecule(a⃗::Atom...) = 
  Molecule([a⃗...])

Molecule(molecule::String) = Molecule(MoleculeName(molecule))
Molecule(::MoleculeName{mn}) where {mn} = error("`molecule` not implemented for $mn")

# parse atomic data from vectors of tuple/vectors
Molecule(atoms::Vector{<:Tuple})  = 
  Molecule([Atom(a[1], a[2:4]) for a in atoms])

Molecule(atoms::Vector{<:Vector}) = 
  Molecule(Tuple.(atoms))

Base.length(molecule::Molecule) = 
  length(atoms(molecule))

Base.getindex(molecule::Molecule, i::Int) = 
  atoms(molecule)[i] 

Base.push!(molecule::Molecule, atom::Atom) = 
  push!(molecule.atoms, atom)

# parse data contained in `mol` to generate
# an input string encoding for HF in Fermi.jl
function parse_molecule(molecule::Molecule)
  molstr = ""
  for a in 1:length(molecule)
    atomname = name(molecule[a])
    atomcoords = Float64.(coordinates(molecule[a]))
    molstr *= atomname 
    for r in atomcoords
      molstr *= " " * string(r)
    end
    molstr *= "\n"
  end
  return molstr
end

function Base.show(io::IO, molecule::Molecule) 
  a⃗ = atoms(molecule)
  println(io, typeof(molecule))
  for (i,a) in enumerate(a⃗)
    print(io, "Atom ",i,": ",name(a))
    x,y,z = coordinates(a)
    print(io, ",   r⃗ = (",x,", ",y,", ",z,")")
    println(io)
  end
end

