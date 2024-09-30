from abc import ABC, abstractmethod
from dataclasses import dataclass

from Bio import PDB
from loguru import logger

from ab_characterisation.developability_tools.tap.definitions import colour_dict


@dataclass
class MetricResult:
    metric_name: str
    calculated_value: float
    flag: str


class BaseMetricCalculator(ABC):
    @abstractmethod
    def __init__(self, quiet: bool = False) -> None:
        self.quiet: bool = quiet
        self.name: str = ""
        self.green_flag_regions: list[tuple[float, float]] = []
        self.amber_flag_regions: list[tuple[float, float]] = []

    def get_flag(self, value: float) -> str:
        """
        Assigns either a green, amber, or red flag to a value depending on the defined regions
        (which were established by calculating the same metrics on structural models of known therapeutics).

        Args:
            value: the value calculated for the query antibody for this metric

        Returns:
            a string representing the flag colour.
        """
        for minval, maxval in self.amber_flag_regions:
            if minval <= value <= maxval:
                return "AMBER"

        for minval, maxval in self.green_flag_regions:
            if minval <= value <= maxval:
                return "GREEN"

        # If the calculated value lies outside the defined green and amber regions, it is assigned a red flag
        return "RED"

    def log_result(self, result: MetricResult) -> None:
        """Logs a message to the terminal summarising the result of the metric calculation."""
        if not self.quiet:
            colour = colour_dict[result.flag]
            logger.opt(colors=True).info(
                f"<b>METRIC {result.metric_name}</b> = <{colour}>{result.calculated_value:.2f} ({result.flag})</{colour}>"
            )

    @abstractmethod
    def calculate(self, annotated_structure: PDB.Structure.Structure) -> MetricResult:
        pass
