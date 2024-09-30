from loguru import logger

from ab_characterisation.developability_tools.sequence_liabilities.scanner_classes import \
    SequenceLiability
from ab_characterisation.developability_tools.utils.outputs import write_file


def display_results(liabilities: list[SequenceLiability]) -> None:
    """
    Nicely prints the identified liabilities to the terminal.

    Args:
        liabilities: the list of identified sequence liabilities.
    """
    color = "red" if liabilities else "green"
    logger.opt(colors=True).info(
        f"\n<b><{color}>{len(liabilities)} liabilities were found</{color}></b>"
    )
    for lia in liabilities:
        logger.opt(colors=True).info(
            f"<b>{lia.liability_type}</b> - residue motif {lia.motif}, position(s) {lia.positions_string}"
        )
    return


def write_liabilities_to_csv(
    liabilities: list[SequenceLiability], filepath: str
) -> None:
    """
    Writes a list of identified sequence liabilities to a file in .csv format.

    Args:
        liabilities: the list of identified sequence liabilities
        filepath: the path to the output file. Can be an S3 path.
    """
    outstr = "Liability,Motif,Positions\n"
    for liability in liabilities:
        outstr += f"{liability.liability_type},{liability.motif},{liability.positions_string}\n"

    write_file(outstr, filepath)

    return
