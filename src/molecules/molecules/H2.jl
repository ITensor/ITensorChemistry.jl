function molecule(::Molecule"H₂")
  return """
  H 0.0 0.0 0.0
  H 0.76 0.0 0.0
  """
end
diis(::Molecule"H₂") = false
oda(::Molecule"H₂") = false
