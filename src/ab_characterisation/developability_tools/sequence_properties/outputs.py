import json

from ab_characterisation.developability_tools.utils.outputs import write_file
from loguru import logger


def write_properties_to_json(property_dict: dict, filepath: str) -> None:
    """
    Writes the calculated property dict to a file in .json format.

    Args:
        property_dict: the dictionary of calculated properties
        filepath: the path to the output file. Can be an S3 path.
    """
    outstr = json.dumps(property_dict)
    write_file(outstr, filepath)

    return
