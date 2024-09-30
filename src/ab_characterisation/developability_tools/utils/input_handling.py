from pathlib import Path
from typing import Optional

from anarci.anarci import number
from Bio import SeqIO


class InputError(Exception):
    pass


def get_numbering(
    sequence: str, expected_type: str
) -> list[tuple[tuple[int, str], str]]:
    """
    Uses ANARCI to number an input sequence.

    Args:
        sequence: the amino acid sequence of the antibody chain
        expected_type: H or L (if the ANARCI annotation does not match this an error will be raised)

    Returns:
        the ANARCI residue numbering, e.g. [((1, ' '), 'E'), ((2, ' '), 'L'), ... ]

    """
    anarci_result: tuple[list[tuple[tuple[int, str], str]], str] = number(sequence)
    numbering, chain_type = anarci_result
    if numbering:
        if chain_type == expected_type:
            return numbering
        raise InputError(
            f"Incorrect chain type: expected {expected_type}, got {chain_type}"
        )
    raise InputError(f"ANARCI failed to number {expected_type} sequence")


def parse_fasta(fasta_file: str) -> dict[str, dict[str, Optional[str]]]:
    if not Path(fasta_file).exists():
        raise InputError(f"Fasta file {fasta_file} does not exist.")

    sequences: dict[str, dict[str, Optional[str]]] = {}
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if '/' not in str(record.seq):
                raise(
                    AssertionError,
                    f"Antibody fasta sequences need to be formatted as HEAVY/LIGHT, an entry in {fasta_file} does not"
                    f"contain /"
                )
            heavy, light = str(record.seq).split("/")
            sequences[record.id] = {
                "H": heavy if heavy not in ["-", ""] else None,
                "L": light if light not in ["-", ""] else None,
            }

    return sequences
