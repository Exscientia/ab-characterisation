from ImmuneBuilder import ABodyBuilder2

from ab_characterisation.developability_tools.tap.main import run_tap as tap
from ab_characterisation.utils.chimerax_utils import ChimeraInput, ChimeraOutput, run_chimerax
from ab_characterisation.utils.data_classes import BiologicsData, RunConfig


def run_abb2(biol_data: BiologicsData, config: RunConfig) -> BiologicsData:
    """
    Run ABodyBuilder2 on input data, save output and output path.
    Args:
        biol_data:
        config:

    Returns:

    """
    predictor = ABodyBuilder2()

    sequences = {"H": biol_data.heavy_sequence, "L": biol_data.light_sequence}

    antibody = predictor.predict(sequences)
    antibody.save(
        str(config.output_directory / "antibody_models" / f"{biol_data.name}_model.pdb")
    )
    biol_data.antibody_structure = (
        config.output_directory / "antibody_models" / f"{biol_data.name}_model.pdb"
    ).resolve()
    return biol_data


def run_tap(biol_data: BiologicsData, config: RunConfig) -> BiologicsData:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    results = tap(biol_data.antibody_structure, outfile=None, quiet=True)
    biol_data.tap_flags = results
    return biol_data


def run_chimerax_superposition(
    biol_data: BiologicsData, config: RunConfig
) -> BiologicsData:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    chimera_input = ChimeraInput(
        name=biol_data.name,
        template=biol_data.target_complex_reference,
        query_ab=biol_data.antibody_structure,
        template_ab_chains=biol_data.target_complex_antibody_chains,
        map_resolution=config.chimera_map_resolution,
        query_ab_chains="HL",
        template_ag_chains=biol_data.target_complex_antigen_chains,
        output_file=str(
            (
                config.output_directory
                / "complex_structures"
                / f"{biol_data.name}_complex.pdb"
            ).resolve()
        ),
    )

    chimera_output = run_chimerax(chimera_input, config)
    if chimera_output.success:
        biol_data.chimerax_complex_structure = chimera_output.output_file
    else:
        biol_data.discarded_by = "ChimeraX failure"

    return biol_data
