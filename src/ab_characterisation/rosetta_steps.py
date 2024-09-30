import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from ab_characterisation.utils.data_classes import BiologicsData, RunConfig


def generic_rosetta_step(
    biol_data: BiologicsData,
    variables: dict[str, str],
    template: str,
    config: RunConfig,
    step_name: str,
    replicates: int = 1,
) -> pd.DataFrame:
    """
    Args:
        biol_data:
        variables:
        template:
        config:
        step_name:
        replicates:

    Returns:

    """
    outputs = []
    for replicate in range(replicates):
        with tempfile.TemporaryDirectory() as temp_dir:
            bash_template_path = (
                Path(__file__).parent / "utils" / "rosetta_templates" / f"{template}.sh"
            )
            xml_template_path = (
                Path(__file__).parent
                / "utils"
                / "rosetta_templates"
                / f"{template}.xml"
            )

            with open(bash_template_path) as inf_sh, open(
                Path(temp_dir) / f"{template}.sh", "w"
            ) as outf_sh:
                for line in inf_sh:
                    for key, value in variables.items():
                        line = line.replace(key, value)
                    outf_sh.write(line)

            shutil.copy(xml_template_path, Path(temp_dir) / f"{template}.xml")

            with open(
                config.output_directory
                / "logs"
                / f"{biol_data.name}_rosetta_{step_name}_{replicate}.log",
                "w",
            ) as outf:
                subprocess.run(
                    ["bash", f"{template}.sh"], cwd=temp_dir, stdout=outf, stderr=outf
                )
            output = pd.read_csv(
                Path(temp_dir) / "score.sc", delim_whitespace=True, skiprows=1
            )
            output["replicate"] = replicate
            outputs.append(output)
    return pd.concat(outputs)


def rosetta_antibody_step(biol_data: BiologicsData, config: RunConfig) -> BiologicsData:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    variables = {
        "<INPUT_FILE>": str(biol_data.antibody_structure),
        "<ROSETTA_BASE_DIR>": config.rosetta_base_directory,
    }
    result_df = generic_rosetta_step(
        biol_data,
        variables,
        "rosetta_metrics_ab_only",
        config,
        step_name="ab_only",
        replicates=config.rosetta_replicates,
    )
    biol_data.rosetta_output_ab_only = result_df
    result_df.to_csv(
        config.output_directory
        / "rosetta_output"
        / f"{biol_data.name}_rosetta_ab_only.csv"
    )
    return biol_data


def rosetta_complex_step(biol_data: BiologicsData, config: RunConfig) -> BiologicsData:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    variables = {
        "<INPUT_FILE>": biol_data.chimerax_complex_structure,
        "<ROSETTA_BASE_DIR>": config.rosetta_base_directory,
    }
    result_df = generic_rosetta_step(
        biol_data,
        variables,
        "rosetta_metrics_complex",
        config,
        step_name="complex",
        replicates=config.rosetta_replicates,
    )
    biol_data.rosetta_output_complex = result_df
    result_df.to_csv(
        config.output_directory
        / "rosetta_output"
        / f"{biol_data.name}_rosetta_complex.csv"
    )
    return biol_data
