from typing import Optional

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class PropertyCalculator:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence
        self.protein_analysis = ProteinAnalysis(self.sequence)

    def calculate_aromaticity(self) -> float:
        return self.protein_analysis.aromaticity()  # type: ignore

    def calculate_charge_at_ph(self, ph: float) -> float:
        return self.protein_analysis.charge_at_pH(ph)  # type: ignore

    def calculate_flexibility(self) -> dict:
        flexibility_scores = self.protein_analysis.flexibility()
        return {
            "residue_scores": flexibility_scores,
            "mean": np.mean(flexibility_scores),
            "stdev": np.std(flexibility_scores),
            "min": min(flexibility_scores),
            "max": max(flexibility_scores),
        }

    def calculate_gravy(self) -> float:
        return self.protein_analysis.gravy()  # type: ignore

    def calculate_instability_index(self) -> float:
        return self.protein_analysis.instability_index()  # type: ignore

    def calculate_isoelectric_point(self) -> float:
        return self.protein_analysis.isoelectric_point()  # type: ignore

    def calculate_properties(self) -> dict:
        """
        Calculates several properties from the sequence of an antibody.
        Returns:
            a dictionary of calculated properties.
        """
        return {
            "sequence": self.sequence,
            "aromaticity": self.calculate_aromaticity(),
            "charge_pH_6": self.calculate_charge_at_ph(6),
            "charge_pH_7.4": self.calculate_charge_at_ph(7.4),
            "flexibility": self.calculate_flexibility(),
            "gravy": self.calculate_gravy(),
            "instability_index": self.calculate_instability_index(),
            "isoelectric_point": self.calculate_isoelectric_point(),
        }
