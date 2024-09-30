import os
from pathlib import Path
from mpi4py import MPI

from ab_characterisation.pipeline_orchestration import pipeline, RunConfig


def test_pipeline():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    input_file = Path(__file__).parent.parent / "data" / "test_pipeline.csv"
    output_dir = Path(__file__).parent.parent / "data" / "ab_characterisation_output"
    rosetta_base_directory = os.environ.get('ROSETTA_BASE')
    config = RunConfig(
        chimera_map_resolution=6,
        input_file=str(input_file),
        output_directory=output_dir,
        rosetta_base_directory=rosetta_base_directory,
    )
    pipeline(config, mpi_rank=rank, mpi_size=size)


if __name__ == '__main__':
    test_pipeline()
