molecule(::MoleculeName"H₂") = 
  molecule([("H", 0.00, 0.00, 0.00),
            ("H", 0.76, 0.00, 0.00)])

molecule(::MoleculeName"N₂") = 
  molecule([("N", 0.00, 0.00,  0.5502960),  
            ("N", 0.00, 0.00, -0.5502960)])

molecule(m::MoleculeName"nitrogen") = molecule(alias(m))
alias(::MoleculeName"nitrogen") = MoleculeName"N₂"()


molecule(::MoleculeName"H₂O") = 
  molecule([("O", 1.2091536548, 1.7664118189, -0.0171613972),  
            ("H", 2.1984800075, 1.7977100627,  0.0121161719),
            ("H", 0.9197881882, 2.4580185570,  0.6297938830)])

molecule(m::MoleculeName"water") = molecule(alias(m))
alias(::MoleculeName"water") = MoleculeName"H₂O"()

molecule(::MoleculeName"BeH₂") = 
  molecule([("Be", 0.00, 0.00,  0.0000),  
            ("H",  0.00, 0.00, -1.3264),
            ("H",  0.00, 0.00, -1.3264)])

molecule(m::MoleculeName"beryllium-dihydride") = molecule(alias(m))
alias(::MoleculeName"beryllium-dihydride") = MoleculeName"BeH₂"()




molecule(::MoleculeName"ammonia") = 
  molecule([("N", 0.0000000,   0.0000000,	 0.1173470), 
            ("H", 0.0000000,   0.9326490,	-0.2738090),
            ("H", 0.8076980,  -0.4663250,	-0.2738090),
            ("H", -0.8076980,	-0.4663250,	-0.2738090)])
molecule(m::MoleculeName"NH₃") = molecule(alias(m))
alias(::MoleculeName"NH₃") = MoleculeName"ammonia"()

molecule(::MoleculeName"formaldehyde") = 
  molecule([("O",	0.00,  0.0000000,  0.6744930),
            ("C",	0.00,  0.0000000, -0.5297240),
            ("H",	0.00,  0.9347280, -1.1087990),
            ("H",	0.00, -0.9347280, -1.1087990)])
molecule(m::MoleculeName"CH₂O") = molecule(alias(m))
alias(::MoleculeName"CH₂O") = MoleculeName"formaldehyde"()


molecule(::MoleculeName"methane") = 
  molecule([("C",	 0.0000000,	 0.0000000,	 0.0000000), 
            ("H",	 0.6268910,	 0.6268910,	 0.6268910),
            ("H",	-0.6268910,	-0.6268910,	 0.6268910),
            ("H",	-0.6268910,	 0.6268910,	-0.6268910),
            ("H",	 0.6268910,	-0.6268910,	-0.6268910)])
molecule(m::MoleculeName"CH₄") = molecule(alias(m))
alias(::MoleculeName"CH₄") = MoleculeName"methane"()

molecule(::MoleculeName"ethanol") = 
  molecule([("C",	 1.1615830,	-0.4067550,	 0.0000000),  
            ("C",	 0.0000000,	 0.5527180,	 0.0000000),
            ("O",	-1.1871140,	-0.2128600,	 0.0000000),
            ("H",	-1.9324340,	 0.3838170,	 0.0000000),
            ("H",	 2.1028600,	 0.1358400,	 0.0000000),
            ("H",	 1.1223470,	-1.0398290,	 0.8811340),
            ("H",	 1.1223470,	-1.0398290,	-0.8811340),
            ("H",	 0.0561470,	 1.1935530,	 0.8808960),
            ("H",	 0.0561470,	 1.1935530,	-0.8808960)])
molecule(m::MoleculeName"C₂H₅OH") = molecule(alias(m))
alias(::MoleculeName"C₂H₅OH") = MoleculeName"ethanol"()


molecule(::MoleculeName"benzene") = 
  molecule([("C",	 0.0000000,	 1.3916730,	0.00), 
            ("C",	 1.2052240,	 0.6958360,	0.00),
            ("C",	 1.2052240,	-0.6958360,	0.00),
            ("C",	 0.0000000,	-1.3916730,	0.00),
            ("C",	-1.2052240,	-0.6958360,	0.00),
            ("C",	-1.2052240,	 0.6958360,	0.00),
            ("H",	 0.0000000,	 2.4695880,	0.00),
            ("H",	 2.1387260,	 1.2347940,	0.00),
            ("H",	 2.1387260,	-1.2347940,	0.00),
            ("H",	 0.0000000,	-2.4695880,	0.00),
            ("H",	-2.1387260,	-1.2347940,	0.00),
            ("H",	-2.1387260,	 1.2347940,	0.00)])
molecule(m::MoleculeName"C₆H₆") = molecule(alias(m))
alias(::MoleculeName"C₆H₆") = MoleculeName"benzene"()

