import typing as t
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ab_characterisation.developability_tools.sequence_liabilities.scanner_classes import \
    SequenceLiability
from ab_characterisation.utils.rosetta_utils import aggregate_rosetta_metrics


@dataclass
class BiologicsData:
    """ """

    heavy_sequence: str
    light_sequence: str
    name: str
    target_complex_reference: str
    target_complex_antigen_chains: str = "A"
    target_complex_antibody_chains: str = "HL"
    antibody_structure: t.Optional[str] = None
    discarded_by: t.Optional[str] = None
    tap_flags: list = field(default_factory=lambda: [])
    sequence_liabilities: list[SequenceLiability] = field(default_factory=lambda: [])
    rosetta_output_ab_only: Optional[pd.DataFrame] = None
    chimerax_complex_structure: t.Optional[str] = None
    rosetta_output_complex: Optional[pd.DataFrame] = None
    rank: Optional[int] = None


@dataclass
class RunConfig:
    """ """

    input_file: str
    output_directory: Path
    rosetta_base_directory: str = None
    chimera_map_resolution: float = 6.0
    dq_sequence_liabilities: list[str] = field(
        default_factory=lambda: ["Unpaired cysteine", "N-linked glycosylation"]
    )
    top_n: int = 100
    rosetta_replicates: int = 1
    exclude_complex_analysis: bool = False

    def __post_init__(self):
        self.output_directory.mkdir(exist_ok=True)
        (self.output_directory / "complex_structures").mkdir(exist_ok=True)
        (self.output_directory / "antibody_models").mkdir(exist_ok=True)
        (self.output_directory / "logs").mkdir(exist_ok=True)
        (self.output_directory / "rosetta_output").mkdir(exist_ok=True)


def save_output(biol_data_ls: list[BiologicsData], config: RunConfig) -> None:
    row_dicts = []
    for biol_data in biol_data_ls:
        row_dict = {}
        for key, value in biol_data.__dict__.items():
            if isinstance(value, str):
                row_dict[key] = value
            elif isinstance(value, int):
                row_dict[key] = value
            elif value is None:
                row_dict[key] = np.nan
            elif key == "tap_flags":
                for tap_metric in value:
                    row_dict[f"TAP-{tap_metric.metric_name}"] = tap_metric.flag
            elif key == "sequence_liabilities":
                seq_liab_str = ""
                for liability in value:
                    seq_liab_str += f"{liability.liability_type}-{liability.motif}-{liability.positions_string}|"
                row_dict[key] = seq_liab_str
            elif key == "rosetta_output_ab_only":
                value = aggregate_rosetta_metrics(value)
                value.columns = ["ab-" + col for col in value.columns]
                for col in value.columns:
                    row_dict[col] = value[col].iloc[0]
            elif key == "rosetta_output_complex":
                value = aggregate_rosetta_metrics(
                    value, metrics=("dG_separated", "total_score")
                )
                value.columns = ["complex-" + col for col in value.columns]
                for col in value.columns:
                    row_dict[col] = value[col].iloc[0]
        row_dicts.append(row_dict)
    pd.DataFrame(row_dicts).to_csv(config.output_directory / "output.csv")
