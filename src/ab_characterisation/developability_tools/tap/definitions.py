colour_dict = {"RED": "red", "AMBER": "fg #ff6100", "GREEN": "fg #14b853"}

imgt_cdr_definitions = {
    1: range(27, 39),
    2: range(56, 66),
    3: range(105, 118),
}

# Two residues on either side of the IMGT CDRs
anchor_residues = [25, 26, 39, 40, 54, 55, 66, 67, 103, 104, 118, 119]


# Salt bridge donor/acceptor atom types
donors = {"LYS": ["NZ"], "ARG": ["NH1", "NH2"]}
acceptors = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}


# Kyte and Doolittle hydrophobicity scale, values normalised to values between 1 and 2
normalised_hydrophobicities = {
    "ILE": 2.0,
    "VAL": 1.9666666666666666,
    "LEU": 1.9222222222222223,
    "PHE": 1.8111111111111111,
    "CYS": 1.7777777777777777,
    "MET": 1.7111111111111112,
    "ALA": 1.7,
    "GLY": 1.4555555555555555,
    "THR": 1.4222222222222223,
    "SER": 1.4111111111111112,
    "TRP": 1.4,
    "TYR": 1.3555555555555556,
    "PRO": 1.3222222222222222,
    "HIS": 1.1444444444444444,
    "GLU": 1.1111111111111112,
    "GLN": 1.1111111111111112,
    "ASP": 1.1111111111111112,
    "ASN": 1.1111111111111112,
    "LYS": 1.0666666666666667,
    "ARG": 1.0,
}

# Charges at pH 7.4
residue_charges = {
    "ALA": 0.0,
    "ARG": 1.0,
    "ASN": 0.0,
    "ASP": -1.0,
    "CYS": 0.0,
    "GLN": 0.0,
    "GLU": -1.0,
    "GLY": 0.0,
    "HIS": 0.1,
    "ILE": 0.0,
    "LEU": 0.0,
    "LYS": 1.0,
    "MET": 0.0,
    "PHE": 0.0,
    "PRO": 0.0,
    "SER": 0.0,
    "THR": 0.0,
    "TRP": 0.0,
    "TYR": 0.0,
    "VAL": 0.0,
}
