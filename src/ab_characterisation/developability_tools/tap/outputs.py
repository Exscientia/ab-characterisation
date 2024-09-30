import tempfile
from pathlib import Path

from ab_characterisation.developability_tools.tap.metrics.base_calculator import MetricResult
from ab_characterisation.developability_tools.utils.outputs import write_file


def write_output_file(results: list[MetricResult], outfile: str) -> None:
    """
    Writes the TAP results to an output file in csv format.

    Args:
        results: the list of metric results
        outfile: the path to where the results should be written.
    """
    outstr = "Metric,Value,Flag\n"
    for res in results:
        outstr += f"{res.metric_name},{res.calculated_value:.2f},{res.flag}\n"

    write_file(outstr, outfile)

    return
