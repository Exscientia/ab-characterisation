from ab_characterisation.developability_tools.sequence_liabilities.scanner_classes import (
    NTerminalGlutamateScanner, RegexScanner)

unpaired_cysteine_scanner = RegexScanner(
    name="Unpaired cysteine",
    description="Checks for C residues in locations other than positions 23 and 104",
    regions=["fv"],
    regex_search_string="C",
    ignored_positions=[(23, " "), (104, " ")],
)

n_linked_glycosylation_scanner = RegexScanner(
    name="N-linked glycosylation",
    description="Checks for an N residue followed by any residue apart from P, followed by S or T",
    regions=["fv"],
    regex_search_string="N[^P][ST]",
)

methionine_oxidation_scanner = RegexScanner(
    name="Methionine oxidation",
    description="Checks for M residues in CDRs or Vernier zones",
    regions=["cdrs", "verniers"],
    regex_search_string="M",
)

tryptophan_oxidation_scanner = RegexScanner(
    name="Tryptophan oxidation",
    description="Checks for W residues in CDRs or Vernier zones",
    regions=["cdrs", "verniers"],
    regex_search_string="W",
)

asparagine_deamidation_scanner = RegexScanner(
    name="Asparagine deamidation",
    description="Checks for residue pairs NG, NS, or NT in CDRs or Vernier zones",
    regions=["cdrs", "verniers"],
    regex_search_string="N[GST]",
)

aspartic_acid_isomeration_scanner = RegexScanner(
    name="Aspartic acid isomeration",
    description="Checks for residue pairs DG, DS, DT, DD, or DH in CDRs or Vernier zones",
    regions=["cdrs", "verniers"],
    regex_search_string="D[GSTDH]",
)

lysine_glycation_scanner = RegexScanner(
    name="Lysine isomeration",
    description="Checks for residue pairs KE, KD, EK, or ED in CDRs or Vernier zones",
    regions=["cdrs", "verniers"],
    regex_search_string="KE|KD|EK|ED",
)

integrin_binding_scanner = RegexScanner(
    name="Integrin binding",
    description="Checks for residue triplets RGD, RYD, or LDV within the Fv",
    regions=["fv"],
    regex_search_string="RGD|RYD|LDV",
)

cd11c_cd18_binding_scanner = RegexScanner(
    name="CD11c/CD18 binding",
    description="Checks for residue triple GPR within the Fv",
    regions=["fv"],
    regex_search_string="GPR",
)

fragmentation_scanner = RegexScanner(
    name="Fragmentation",
    description="Checks for residue pair DP in the CDRs or Vernier zones",
    regions=["cdrs", "verniers"],
    regex_search_string="DP",
)

n_terminal_glutamate_scanner = NTerminalGlutamateScanner(
    name="N-terminal glutamate",
    description="Checks for glutamate residues at the N-termini of both heavy and light chains",
)
