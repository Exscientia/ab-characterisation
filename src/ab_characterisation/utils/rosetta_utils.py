import pandas as pd


def aggregate_rosetta_metrics(
    metric_df: pd.DataFrame, metrics: list[str] = ["dG_separated"]
) -> pd.DataFrame:
    """
    Args:
        metric_df:
        metrics:

    Returns:

    """
    metric_df = metric_df.select_dtypes("number")
    if len(metrics) == 1:
        idx = list(metric_df.sort_values(by=metrics[0], ascending=True)[:3].index)

    else:
        idx = []
        for metric in metrics:
            idx += list(metric_df.sort_values(by=metric, ascending=True)[:2].index)
        idx = set(idx)
    metric_df = metric_df.loc[list(idx)].mean().to_frame().T
    return metric_df
