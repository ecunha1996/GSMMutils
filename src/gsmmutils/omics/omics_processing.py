from typing import Union

import numpy as np
import pandas
import pandas as pd
from numpy import log


def global_thresholding(sample_series, gtu, maxexp):
    return (sample_series / gtu).apply(log).clip(-maxexp, maxexp).to_dict()


def local1_thresholding(sample_series, gtu, lt, maxexp):
    activity = (sample_series / gtu).apply(log).clip(-maxexp, maxexp)
    gt_active = activity >= 0
    activity[gt_active] = (log(sample_series[gt_active] / lt[gt_active]) * maxexp).clip(-maxexp, maxexp)
    return activity.to_dict()


def local2_thresholding(sample_series, gtu, gtl, lt, maxexp):
    upp_activity = (1 + (sample_series / gtu).apply(log)).clip(-maxexp, 1 + maxexp)
    gtu_inactive = upp_activity < 1
    low_activity = (sample_series / gtl).apply(log).clip(-maxexp, maxexp)
    gtl_maybes, gtl_lows = (low_activity >= 0) & gtu_inactive, low_activity < 0
    upp_activity[gtl_lows] = low_activity[gtl_lows]
    activity_maybe = (sample_series[gtl_maybes] / lt[gtl_maybes]). \
        apply(log).clip(-maxexp, maxexp)
    upp_activity[gtl_maybes] = activity_maybe.clip(-1, 1)
    return upp_activity.to_dict()


def threshold_strategy(sample_series, gtu, gtl, lt, maxexp, thresholding_strategy) -> dict:
    """
    Thresholding strategy for the omics data. Processes a single sample at the time.

    Parameters
    ----------
    sample_series:
        Omics data from a specific sample.
    thresholding_strategy:
        Thresholding strategy to be used.
    gtu:
        Global threshold upper value.
    gtl:
        Global threshold lower value.
    lt:
        Local threshold value.
    maxexp
        Maximum expression value of the dataset.

    Returns
    -------
    filtered_sample: dict
        Filtered omics data from a specific sample.

    """
    if thresholding_strategy == 'Global':
        return global_thresholding(sample_series, gtu, maxexp)

    elif thresholding_strategy == 'Local1':
        return local1_thresholding(sample_series, gtu, lt, maxexp)

    elif thresholding_strategy == 'Local2':
        return local2_thresholding(sample_series, gtu, gtl, lt, maxexp)

    else:
        raise ValueError('Thresholding strategy not recognized')


def thresholding_filter(omics_dataframe: pandas.DataFrame, thresholding_strategy: Union[int, str],
                        global_threshold_upper: int, global_threshold_lower: int,
                        local_threshold: int) -> pandas.DataFrame:
    """
    Thresholding filter for the omics data.

    Parameters
    ----------
    omics_dataframe: pandas.DataFrame
        Omics dataframe to be filtered.
    thresholding_strategy: str
        Thresholding strategy to be used.
        Accepted values: 'global', 'local1', 'local2'
    global_threshold_upper: int
        Position of the Global threshold value on the quantile list.
    global_threshold_lower: int
        Position of the Local threshold 1 value on the quantile list.
    local_threshold: int
        Position of the Local threshold 2 value on the quantile list.

    Returns
    -------
    filtered_dataset: pandas.DataFrame
        Filtered omics dataframe.

    """
    max_expression = np.log(omics_dataframe.max().max())

    qvalues = [0.1, 0.25, 0.5, 0.75, 0.9]
    quantiles = omics_dataframe.quantile(qvalues)
    global_thresholds = quantiles.T.apply(lambda x: x.mean())

    filter_results = {}

    for sample_id in omics_dataframe.index:
        sample = omics_dataframe.loc[sample_id, :]

        name = '_'.join(map(str, [thresholding_strategy, global_threshold_upper,
                                  global_threshold_lower, local_threshold]))

        gtl = global_thresholds.iloc[global_threshold_upper]
        gtu = global_thresholds.iloc[global_threshold_lower]
        lti = quantiles.iloc[local_threshold, ]

        filter_results[sample.name + '_' + name] = threshold_strategy(sample, gtl, gtu, lti, max_expression,
                                                                      thresholding_strategy)

    filtered_dataset = pd.DataFrame(filter_results).T

    return filtered_dataset
