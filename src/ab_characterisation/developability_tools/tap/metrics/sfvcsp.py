from Bio import PDB

from ab_characterisation.developability_tools.tap.metrics.base_calculator import (
    BaseMetricCalculator, MetricResult)


class SFvCSPCalculator(BaseMetricCalculator):
    def __init__(self, quiet: bool = False) -> None:
        self.quiet = quiet
        self.name = "SFvCSP"
        self.green_flag_regions = [(-4.20, 100000)]
        self.amber_flag_regions = [(-20.50, -4.20)]

    def calculate(self, annotated_structure: PDB.Structure.Structure) -> MetricResult:
        """
        Calculates the 'structural Fv charge symmetry parameter' (SFvCSP) and assigns a flag colour.
        Considers surface residues only.
        The input structure must have been annotated using the
        ab_characterisation.developability_tools.tap.structure_annotation module.
        """
        h_charge = sum(
            res.charge
            for res in annotated_structure[0]["H"].get_residues()
            if res.is_surface
        )
        l_charge = sum(
            res.charge
            for res in annotated_structure[0]["L"].get_residues()
            if res.is_surface
        )
        sfvcsp = h_charge * l_charge
        flag = self.get_flag(sfvcsp)

        result = MetricResult(
            metric_name=self.name,
            calculated_value=sfvcsp,
            flag=flag,
        )

        self.log_result(result)
        return result
