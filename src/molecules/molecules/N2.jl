function molecule(::Molecule"N₂")
  return """
  N	0.0000000	0.0000000	0.5502960
  N	0.0000000	0.0000000	-0.5502960
  """
end
molecule(m::Molecule"nitrogen") = molecule(alias(m))
alias(::Molecule"nitrogen") = Molecule"N₂"()
