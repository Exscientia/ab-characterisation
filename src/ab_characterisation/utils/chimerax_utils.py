import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from ImmuneBuilder.refine import refine

from ab_characterisation.utils.data_classes import RunConfig


@dataclass
class ChimeraInput:
    name: str
    template: str
    query_ab: str
    template_ab_chains: str
    map_resolution: float
    query_ab_chains: str
    template_ag_chains: str
    output_file: str


@dataclass
class ChimeraOutput:
    success: bool
    output_file: str


def write_script(script_name: str, payload: ChimeraInput) -> None:
    """

    Args:
        script_name:
        payload:

    Returns:

    """
    with open(script_name, "w") as outf:
        outf.write("from chimerax.core.commands import run\n")
        outf.write(f"run(session, 'open {payload.template}')\n")
        outf.write(
            f"run(session, 'molmap /{','.join(list(payload.template_ab_chains))} {payload.map_resolution}')\n"
        )
        outf.write(f"run(session, 'open {payload.query_ab}')\n")
        outf.write("run(session, 'fitmap #3 inMap #2 search 10')\n")
        outf.write(
            f"""run(session, "select #3/{','.join(list(payload.query_ab_chains))}#1/{','.join(list(payload.template_ag_chains))}")\n"""
        )
        outf.write(
            f"""run(session, "save {payload.output_file} format pdb selectedOnly true")\n"""
        )
        outf.write("""run(session, "exit")\n""")


def run_chimerax(payload: ChimeraInput, config: RunConfig) -> ChimeraOutput:
    """
    Use Chimerax to create complex pdb file of the query AB and the target antigen, using the template context to guide
    the complex generation.
    Args:
        payload:

    Returns:

    """
    with tempfile.NamedTemporaryFile(suffix=".py") as temp_f:
        script_name = temp_f.name
        write_script(payload=payload, script_name=script_name)

        cmd = ["ChimeraX", "--script", script_name, "--nogui"]
        with open(
            config.output_directory / "logs" / f"{payload.name}_chimera.log", "w"
        ) as outf:
            subprocess.run(cmd, check=True, stderr=outf, stdout=outf)

    if not Path(payload.output_file).exists():
        output = ChimeraOutput(output_file=payload.output_file, success=False)
        return output

    with open(payload.output_file) as inf:
        lines = inf.readlines()
    with open(payload.output_file, "w") as outf:
        for line in lines:
            if line.startswith("ATOM"):
                outf.write(line)

    # refinement
    refined_output = payload.output_file.replace(".pdb", "_refined.pdb")
    success = refine(input_file=payload.output_file, output_file=refined_output)
    if not success:
        output = ChimeraOutput(output_file=refined_output, success=success)
        return output

    output = ChimeraOutput(
        output_file=payload.output_file, success=Path(payload.output_file).exists()
    )
    return output
