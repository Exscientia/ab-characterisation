from typing import Optional

from loguru import logger
import numpy as np
import pandas as pd
from numpy import typing as npt
from scipy.stats import multivariate_normal

from ab_characterisation.utils.data_classes import BiologicsData, RunConfig
from ab_characterisation.utils.rosetta_utils import aggregate_rosetta_metrics


def sequence_liability_filter(biol_data: BiologicsData, config: RunConfig) -> bool:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    for liability in biol_data.sequence_liabilities:
        if liability.liability_type in config.dq_sequence_liabilities:
            return True
    return False


def tap_filter(biol_data: BiologicsData, config: RunConfig) -> bool:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    for tap_metric in biol_data.tap_flags:
        if tap_metric.flag == "RED":
            return True
    return False


def rosetta_antibody_filter(biol_data: BiologicsData, config: RunConfig) -> bool:
    """

    Args:
        biol_data:
        config:

    Returns:

    """
    rosetta_antibody_data = aggregate_rosetta_metrics(biol_data.rosetta_output_ab_only)
    if rosetta_antibody_data.dG_separated.iloc[0] < 5:
        return False
    return True


def find_top_n(
    biol_data_ls: list[BiologicsData], config: RunConfig
) -> list[BiologicsData]:
    pd_row_ls = []
    out_ls = []
    for biol_data in biol_data_ls:
        if biol_data.discarded_by is None:
            rosetta_complex_scores = aggregate_rosetta_metrics(
                biol_data.rosetta_output_complex, metrics=["dG_separated", "total_score"]
            ).iloc[0]
            pd_row_ls.append(rosetta_complex_scores)
    metric_df = pd.concat(pd_row_ls)
    top_indices = find_top_candidates(
        metric_df.total_score,
        metric_df.dG_separated,
        config.top_n,
        scale_factor=1,
        fit_without_outliers=True,
    )
    valid_candidate_idx = -1
    for biol_data in biol_data_ls:
        if biol_data.discarded_by is None:
            valid_candidate_idx += 1
            if valid_candidate_idx in top_indices:
                biol_data.rank = list(top_indices).index(valid_candidate_idx)
            else:
                biol_data.discarded_by = "Not in top N"
        out_ls.append(biol_data)
    return out_ls


def find_top_candidates(
    total_score: npt.ArrayLike,
    dG_separated: npt.ArrayLike,
    n: int,
    scale_factor: float = 1,
    total_score_max: Optional[float] = None,
    dG_separated_max: Optional[float] = None,
    fit_without_outliers: bool = True,
) -> np.ndarray:
    """
    Fit multivariate gaussian (centered on median rather than mean) and then select best points
    according to lowest probability of being drawn subject to bounds.
    Args:
        total_score: Total score for candidates to select
        dG_separated: dG_seperated of candidates to select
        n: N candidates to select
        scale_factor: Scale total score of data points by this factor AFTER fitting multivariate
        total_score_max: Maximum total score of selected candidates, if not specified use median
        dG_separated_max: Maximum dG_separated of selected candidates, if not specified use medians
        fit_without_outliers: Ignore points that are 1.5 IQR above/below the upper/lower quartile
                              when fitting the gaussian.

    Returns:

    """
    data = np.stack([np.array(total_score), np.array(dG_separated)], axis=1)

    if fit_without_outliers:
        # Define outliers (don't fit gaussian on these)
        ts_q1 = np.quantile(total_score, 0.25)
        ts_q3 = np.quantile(total_score, 0.75)
        ts_IQR = ts_q3 - ts_q1
        ts_lower_bound = ts_q1 - ts_IQR * 1.5
        ts_upper_bound = ts_q3 + ts_IQR * 1.5

        dGs_q1 = np.quantile(dG_separated, 0.25)
        dGs_q3 = np.quantile(dG_separated, 0.75)
        dGs_IQR = dGs_q3 - dGs_q1
        dGs_lower_bound = dGs_q1 - dGs_IQR * 1.5
        dGs_upper_bound = dGs_q3 + dGs_IQR * 1.5

        outlier_mask = (
            (ts_lower_bound < data[:, 0])
            & (data[:, 0] < ts_upper_bound)
            & (dGs_lower_bound < data[:, 1])
            & (data[:, 1] < dGs_upper_bound)
        )

        median = np.median(data[outlier_mask], axis=0)
        cov = np.cov(data[outlier_mask], rowvar=0)
    else:
        median = np.median(data, axis=0)
        cov = np.cov(data, rowvar=0)
    multivar_f = multivariate_normal(mean=median, cov=cov, allow_singular=True)
    xmax = total_score_max if total_score_max is not None else median[0]
    ymax = dG_separated_max if dG_separated_max is not None else median[1]
    centroid = np.array([xmax, ymax])
    mask = np.all((data < centroid), axis=1)
    data[:, 0] = scale_factor * (data[:, 0] - median[0]) + median[0]
    top_idx = np.argsort(multivar_f.pdf(data[mask]))[:n]
    indices = np.where(mask == True)[0]
    return indices[top_idx]
