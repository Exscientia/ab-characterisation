from ab_characterisation.developability_tools.sequence_liabilities.main import scan_single
from ab_characterisation.utils.data_classes import BiologicsData, RunConfig


def sequence_liability_check(
    input_data: BiologicsData, config: RunConfig
) -> BiologicsData:
    """

    Args:
        input_data:

    Returns:

    """
    liabilities = scan_single(
        input_data.heavy_sequence, input_data.light_sequence, quiet=True
    )
    input_data.sequence_liabilities = liabilities
    return input_data
