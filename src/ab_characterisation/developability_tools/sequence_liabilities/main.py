import sys
from typing import Optional

from loguru import logger

from ab_characterisation.developability_tools.utils.input_handling import get_numbering
from ab_characterisation.developability_tools.sequence_liabilities.scanner_classes import SequenceLiability
from ab_characterisation.developability_tools.sequence_liabilities.scanners import (asparagine_deamidation_scanner,
                                                                    aspartic_acid_isomeration_scanner,
                                                                    cd11c_cd18_binding_scanner, fragmentation_scanner,
                                                                    integrin_binding_scanner, lysine_glycation_scanner,
                                                                    methionine_oxidation_scanner,
                                                                    n_linked_glycosylation_scanner,
                                                                    n_terminal_glutamate_scanner,
                                                                    tryptophan_oxidation_scanner, unpaired_cysteine_scanner)

logger.remove()
logger.add(sys.stderr, format="{message}")


scanner_list = [
    unpaired_cysteine_scanner,
    n_linked_glycosylation_scanner,
    methionine_oxidation_scanner,
    tryptophan_oxidation_scanner,
    asparagine_deamidation_scanner,
    aspartic_acid_isomeration_scanner,
    lysine_glycation_scanner,
    integrin_binding_scanner,
    cd11c_cd18_binding_scanner,
    fragmentation_scanner,
    n_terminal_glutamate_scanner,
]


def scan_single(
    heavy_sequence: Optional[str], light_sequence: Optional[str], quiet: bool = False
) -> list[SequenceLiability]:
    """
    Scans the sequence of an antibody for potential liabilities.

    Args:
        heavy_sequence: the amino acid sequence of the antibody heavy chain
        light_sequence: the amino acid sequence of the antibody light chain

    Returns:
        a list of identified sequence liabilities.
    """

    numbering_dict = {}
    if heavy_sequence:
        numbering_dict["H"] = get_numbering(heavy_sequence, "H")
    if light_sequence:
        numbering_dict["L"] = get_numbering(light_sequence, "L")

    liabilities = []
    for scanner in scanner_list:
        liabilities += scanner.scan(numbering_dict, quiet=quiet)

    return liabilities
