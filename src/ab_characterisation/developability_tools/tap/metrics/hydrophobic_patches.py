from Bio import PDB

from ab_characterisation.developability_tools.tap.metrics.base_calculator import (
    BaseMetricCalculator, MetricResult)


class HydrophobicPatchScoreCalculator(BaseMetricCalculator):
    def __init__(self, quiet: bool = False) -> None:
        self.quiet = quiet
        self.name = "Hydrophobic Patch Score"
        self.green_flag_regions = [(137.61, 200.71)]
        self.amber_flag_regions = [(106.44, 137.61), (200.71, 225.85)]

    def calculate(self, annotated_structure: PDB.Structure.Structure) -> MetricResult:
        """
        Calculates the 'patches of surface hydrophobicity' score (PSH) and assigns a flag colour.
        Considers residues that are in the CDR vicinity only.
        The input structure must have been annotated using the
        ab_characterisation.developability_tools.tap.structure_annotation module.
        """
        cdr_vicinity = [
            res for res in annotated_structure[0].get_residues() if res.in_cdr_vicinity
        ]

        score = 0
        for res1 in cdr_vicinity:
            for res2 in cdr_vicinity:
                if res1 == res2:
                    continue

                distance = res1.neighbours.get(res2.get_full_id(), None)
                if not distance:
                    continue

                score += (res1.hydrophobicity * res2.hydrophobicity) / distance**2

        flag = self.get_flag(score)

        result = MetricResult(
            metric_name=self.name,
            calculated_value=score,
            flag=flag,
        )

        self.log_result(result)
        return result
