from Bio import PDB

from ab_characterisation.developability_tools.tap.metrics.base_calculator import (
    BaseMetricCalculator, MetricResult)


class TotalCDRLengthCalculator(BaseMetricCalculator):
    def __init__(self, quiet: bool = False) -> None:
        self.quiet = quiet
        self.name = "Total IMGT CDR Length"
        self.green_flag_regions = [(43, 55)]
        self.amber_flag_regions = [(37, 43), (55, 63)]

    def calculate(self, annotated_structure: PDB.Structure.Structure) -> MetricResult:
        """
        Calculates the total number of CDR residues and assigns a flag colour.
        Uses the IMGT CDR definition.
        The input structure must have been annotated using the
        ab_characterisation.developability_tools.tap.structure_annotation module.
        """
        cdr_residues = [
            res for res in annotated_structure[0].get_residues() if res.is_cdr
        ]
        total_cdr_length = len(cdr_residues)
        flag = self.get_flag(total_cdr_length)

        result = MetricResult(
            metric_name=self.name,
            calculated_value=total_cdr_length,
            flag=flag,
        )

        self.log_result(result)
        return result
