import sys
import typing as t

import pandas as pd
from loguru import logger
from mpi4py import MPI
from ab_characterisation.utils.data_classes import BiologicsData, RunConfig, save_output

from ab_characterisation.filter_steps import (
    find_top_n, rosetta_antibody_filter, sequence_liability_filter, tap_filter
)
from ab_characterisation.rosetta_steps import rosetta_antibody_step, rosetta_complex_step
from ab_characterisation.sequence_steps import sequence_liability_check
from ab_characterisation.structure_steps import run_abb2, run_chimerax_superposition, run_tap


def get_objects(config: RunConfig) -> list[BiologicsData]:
    """

    Args:
        config:

    Returns:

    """
    data_objects = []
    df = pd.read_csv(config.input_file)
    for idx, row in df.iterrows():
        data_objects.append(
            BiologicsData(
                heavy_sequence=row.heavy_sequence,
                light_sequence=row.light_sequence,
                name=row.sequence_name,
                target_complex_reference=row.reference_complex,
            )
        )
    return data_objects


def filtering_step(
    input_data: list[BiologicsData],
    step_name: str,
    criterion_function: t.Callable,
    config: RunConfig,
) -> list[BiologicsData]:
    """
    General framework for a step that performs filtering of the input data, labelling datapoints as discarded if they
    fail to pass a filter criterion.
    Args:
        input_data:
        step_name:
        criterion_function: Function mapping BiologicsData -> bool
        config:

    Returns:
        list of BiologicsData objects
    """
    output_data: list[BiologicsData] = []
    filter_count = 0

    for biol_data in input_data:
        if biol_data.discarded_by is None:
            filtered = criterion_function(biol_data, config)
            if filtered:
                biol_data.discarded_by = step_name
                filter_count += 1
        output_data.append(biol_data)

    logger.info(f"{filter_count} datapoints discarded during step {step_name}.")
    return output_data


def computation_step(
    input_data: list[BiologicsData], computation_function: t.Callable, config: RunConfig
) -> list[BiologicsData]:
    """
    General framework for a step that performs computation on the input data, manipulating one or more of the dataclass
    fields.

    Args:
        input_data:
        computation_function: Function mapping BiologicsData -> BiologicsData, modifying the dataclass fields with the
            results of the computation
        config:

    Returns:
        list of BiologicsData objects
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Calculate the chunk size for each process
    chunk_size = len(input_data) // size
    remainder = len(input_data) % size

    # Calculate the range for the current process
    local_start = rank * chunk_size + min(rank, remainder)
    local_end = local_start + chunk_size + (1 if rank < remainder else 0)

    # Perform the local computation
    local_results: list[BiologicsData] = []
    for biol_data in input_data[local_start:local_end]:
        if biol_data.discarded_by is None:
            biol_data = computation_function(biol_data, config)
        local_results.append(biol_data)

    # Gather the local results at the root process
    all_results = comm.gather(local_results, root=0)

    # Combine the results into a single list
    output_data: list[BiologicsData] = []
    if rank == 0:
        for result_list in all_results:
            output_data.extend(result_list)
    output_data = comm.bcast(output_data, root=0)
    return output_data


def pipeline(config: RunConfig, mpi_rank: int, mpi_size: int) -> None:
    """

    Args:
        config:

    Returns:

    """
    logger.remove()
    if mpi_rank == 0:
        logger.add(
            sys.stdout,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
            level="INFO",
        )
    else:
        logger.add(
            sys.stdout,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
            level="WARNING",
        )
    biologics_objects = get_objects(config)

    logger.info("Identifying sequence liabilities")
    biologics_objects = computation_step(
        biologics_objects, sequence_liability_check, config
    )
    logger.info("Filtering by sequence liabilities")
    biologics_objects = filtering_step(
        biologics_objects,
        step_name="liabilities",
        criterion_function=sequence_liability_filter,
        config=config,
    )

    logger.info("Running ABB2")
    biologics_objects = computation_step(biologics_objects, run_abb2, config)
    logger.info("Running TAP")
    biologics_objects = computation_step(biologics_objects, run_tap, config)
    logger.info("Filtering TAP")
    biologics_objects = filtering_step(biologics_objects, "tap", tap_filter, config)
    logger.info("Running antibody-only Rosetta analysis")
    biologics_objects = computation_step(
        biologics_objects, rosetta_antibody_step, config
    )
    logger.info("Running filtering based on antibody-only Rosetta analysis")
    biologics_objects = filtering_step(
        biologics_objects, "rosetta_antibody", rosetta_antibody_filter, config
    )
    if not config.exclude_complex_analysis:
        logger.info("Running ChimeraX complex generation")
        biologics_objects = computation_step(
            biologics_objects, run_chimerax_superposition, config
        )
        logger.info("Running Rosetta complex analysis")
        biologics_objects = computation_step(
            biologics_objects, rosetta_complex_step, config
        )

    if mpi_rank == 0:
        logger.info("Identifying top N candidates")
        biologics_objects = find_top_n(biologics_objects, config)
        save_output(biol_data_ls=biologics_objects, config=config)
