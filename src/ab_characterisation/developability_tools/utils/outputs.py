import tempfile
from pathlib import Path


def write_file(contents: str, filepath: str) -> None:
    """Writes an output file to the given location with the given contents."""
    outpath = Path(filepath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w") as openf:
        openf.write(contents)
