from mpi4py import MPI
from pathlib import Path
import typer

from ab_characterisation.pipeline_orchestration import RunConfig, pipeline

app = typer.Typer(
    name="ab-characterisation-pipeline",
    add_completion=False,
)


@app.command()
def run_pipeline(
    input_file: str = typer.Option(..., help='Input .csv file, containing sequence_name, heavy_sequence, light_sequence '
                                             'and reference_complex columns.'),
    chimera_resolution: float = typer.Option(6.0, help='Resolution of the map used for alignment within ChimeraX.'),
    output_dir: str = typer.Option("./ab_characterisation_output", help='Directory to which output files are written.'),
    rosetta_replicates: int = typer.Option(1, help='How many replicates to run for Rosetta characterisation steps.'),
    rosetta_base_dir: str = typer.Option(..., help='Base directory for the Roestta software suite, e.g. '
                                                   '/path/to/rosetta/rosetta.binary.linux.release-315'),
    top_n: int = typer.Option(10, help='Top N candidate antibodies to provide from the provided .csv file of antibodies'),
    no_complex_analysis: bool = typer.Option(False, help='If provided, the pipeline does not perform antibody-antigen '
                                                         'complex generation and analysis.')
):
    output_dir = Path(output_dir)
    config = RunConfig(
        chimera_map_resolution=chimera_resolution,
        input_file=input_file,
        output_directory=output_dir,
        rosetta_base_directory=rosetta_base_dir,
        top_n=top_n,
        rosetta_replicates=rosetta_replicates,
        exclude_complex_analysis=no_complex_analysis,
    )
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    pipeline(config, mpi_rank=rank, mpi_size=size)


if __name__ == "__main__":
    app()

