import sys
from pathlib import Path
from typing import Optional

from ab_characterisation.developability_tools.sequence_properties.calculations import \
    PropertyCalculator
from ab_characterisation.developability_tools.sequence_properties.outputs import \
    write_properties_to_json
from ab_characterisation.developability_tools.utils.input_handling import parse_fasta
from loguru import logger

logger.remove()
logger.add(sys.stderr, format="{message}")


def calculate_properties(
    heavy_sequence: Optional[str], light_sequence: Optional[str]
) -> dict[str, dict]:
    results = {}
    for chain, sequence in {"heavy": heavy_sequence, "light": light_sequence}.items():
        if sequence:
            results[chain] = PropertyCalculator(sequence).calculate_properties()

    return results


def property_calculator(
    heavy_sequence: Optional[str],
    light_sequence: Optional[str],
    outfile: Optional[str] = None,
) -> dict[str, dict]:
    """
    Main function to calculate sequence-based properties for a single antibody.
    Writes the results to a file in .csv format.

    Args:
        heavy_sequence: the amino acid sequence of the antibody heavy chain
        light_sequence: the amino acid sequence of the antibody light chain
        outfile: the output path where results should be written. If a path is not given, results will only be printed
            to the terminal.
    """
    property_dict = calculate_properties(heavy_sequence, light_sequence)

    if outfile:
        write_properties_to_json(property_dict, outfile)

    return property_dict


def property_calculator_fasta(
    fasta_file: str,
    outdir: Optional[str],
    quiet: bool = False,
) -> dict[str, dict[str, dict]]:
    """
    Function to scan a set of antibody sequences in a fasta file for liabilities.
    Writes the results to a series of files (one per antibody) in .csv format.

    Args:
        fasta_file: the amino acid sequences of the antibodies to be scanned in fasta format. Each antibody should be a
            separate fasta entry, with the heavy and light chains being separated by a forward slash. E.g.:

                    >antibody1
                    HEAVYSEQUENCE/LIGHTSEQUENCE
                    >antibody2
                    HEAVYSEQUENCE/LIGHTSEQUENCE
                    >nanobody1
                    HEAVYSEQUENCE/-
                    ...

        outdir: Path to the directory where results should be written. Individual files will be named according to the
            IDs in the fasta file.
    """
    if outdir:
        dirpath = Path(outdir)
        dirpath.mkdir(parents=True, exist_ok=True)
    else:
        dirpath = Path(".")

    sequences = parse_fasta(fasta_file)
    results = {}
    for antibody_id, seqs in sequences.items():
        if not quiet:
            logger.info(f"Calculating properties for {antibody_id}")

        property_dict = calculate_properties(seqs["H"], seqs["L"])
        filepath = dirpath / f"{antibody_id}_properties.json"
        write_properties_to_json(property_dict, str(filepath))
        results[antibody_id] = property_dict

    return results
