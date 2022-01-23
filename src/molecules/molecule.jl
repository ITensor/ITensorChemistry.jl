struct Molecule{M} end

Molecule(molecule::Symbol) = Molecule{molecule}()
Molecule(molecule::String) = Molecule(Symbol(molecule))
macro Molecule_str(s)
    :(Molecule{$(Expr(:quote, Symbol(s)))})
end

# Interface
molecule(::Molecule{M}) where {M} = error("`molecule` not implemented for $M")
diis(::Molecule) = true
oda(::Molecule) = true
