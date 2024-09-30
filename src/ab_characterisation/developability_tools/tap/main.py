import sys
from typing import Optional

from loguru import logger

from ab_characterisation.developability_tools.tap.definitions import colour_dict
from ab_characterisation.developability_tools.tap.metrics import (
    HydrophobicPatchScoreCalculator, NegativePatchScoreCalculator,
    PositivePatchScoreCalculator, SFvCSPCalculator, TotalCDRLengthCalculator)
from ab_characterisation.developability_tools.tap.metrics.base_calculator import MetricResult
from ab_characterisation.developability_tools.tap.outputs import write_output_file
from ab_characterisation.developability_tools.tap.structure_annotation import \
    StructureAnnotator

logger.remove()
logger.add(sys.stderr, format="{message}")


def run_tap(
    modelfile: str,
    outfile: Optional[str],
    quiet: bool = False,
) -> list[MetricResult]:
    """
    Main function to calculate TAP metrics for a pre-generated ABodyBuilder2 model.
    Writes the results to a file in .csv format.

    Args:
        modelfile: the path to the input model .pdb file. Can be an S3 path - in this case the file will be temporarily
            downloaded before TAP is run.
            This should be a model created by ABodyBuilder2, and should be already IMGT numbered.
        outfile: the output path where results should be written.
        quiet: suppresses all log messages if set to True.
    """

    structure = StructureAnnotator().load_and_annotate_structure(modelfile)

    # Calculate the 5 metrics
    results = []
    for calculator in [
        HydrophobicPatchScoreCalculator,
        NegativePatchScoreCalculator,
        PositivePatchScoreCalculator,
        SFvCSPCalculator,
        TotalCDRLengthCalculator,
    ]:
        results.append(calculator(quiet=quiet).calculate(structure))  # type: ignore

    if outfile:
        write_output_file(results, outfile)

    return results


def list_metrics() -> list[dict]:
    """
    Returns (and prints) a list of the metrics and their green/amber region definitions.
    """
    metrics = []
    for calculator in [
        HydrophobicPatchScoreCalculator,
        NegativePatchScoreCalculator,
        PositivePatchScoreCalculator,
        SFvCSPCalculator,
        TotalCDRLengthCalculator,
    ]:
        metric = calculator()  # type: ignore
        green_str = "; ".join(
            [
                str(region[0]) + " to " + str(region[1])
                for region in metric.green_flag_regions
            ]
        )
        amber_str = "; ".join(
            [
                str(region[0]) + " to " + str(region[1])
                for region in metric.amber_flag_regions
            ]
        )
        logger.opt(colors=True).info(f"<b>TAP METRIC {metric.name}:</b>")
        logger.opt(colors=True).info(
            f"<{colour_dict['GREEN']}>GREEN</{colour_dict['GREEN']}> region: {green_str}"
        )
        logger.opt(colors=True).info(
            f"<{colour_dict['AMBER']}>AMBER</{colour_dict['AMBER']}> region: {amber_str}\n"
        )

        metrics.append(
            {
                "name": metric.name,
                "green_flag_regions": metric.green_flag_regions,
                "amber_flag_regions": metric.amber_flag_regions,
            }
        )

    return metrics
