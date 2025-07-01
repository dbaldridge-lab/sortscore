"""
Score calculation functions for variant activity analysis.

This module provides functions to calculate activity scores from count data, supporting different averaging methods and filtering criteria.

Examples
--------
>>> from sortscore.analysis.score import calculate_activity_scores
"""
import logging
import pandas as pd
from typing import List, Optional

def merge_count_dfs(count_dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge multiple count DataFrames on the 'variant' column.
    """
    merged = count_dfs[0]
    for df in count_dfs[1:]:
        merged = pd.merge(merged, df, on='variant', how='outer', suffixes=(False, False))
    merged = merged.fillna(0)
    return merged

def normalize_counts(merged: pd.DataFrame, read_count: Optional[List[int]]) -> pd.DataFrame:
    """
    Normalize count columns by total reads per sample.
    """
    if read_count is not None:
        for i, rc in enumerate(read_count):
            col = merged.columns[i+1]  # skip 'variant' column
            merged[col] = merged[col] / rc * 1e6
    return merged

def filter_variants(
    merged: pd.DataFrame,
    min_reads: int,
    bins_required: int,
    reps_required: int
) -> pd.DataFrame:
    """
    Filter variants by minimum reads, bins, and replicates required.
    """
    count_cols = merged.columns[1:]
    # Filter by min_reads in all columns
    if min_reads > 0:
        merged = merged[(merged[count_cols] >= min_reads).all(axis=1)]
    # Filter by bins_required (number of columns with nonzero counts)
    merged['bins_present'] = (merged[count_cols] > 0).sum(axis=1)
    merged = merged[merged['bins_present'] >= bins_required]
    # Filter by reps_required (number of columns with nonzero counts)
    merged['reps_present'] = (merged[count_cols] > 0).sum(axis=1)
    merged = merged[merged['reps_present'] >= reps_required]
    merged = merged.drop(columns=['bins_present', 'reps_present'])
    return merged

def simple_average(merged: pd.DataFrame) -> pd.Series:
    """
    Calculate the simple average across all count columns.
    """
    return merged.iloc[:, 1:].mean(axis=1)

def rep_weighted_average(merged: pd.DataFrame, read_count: Optional[List[int]]) -> pd.Series:
    """
    Calculate replicate-weighted average across count columns.
    """
    count_cols = merged.columns[1:]
    if read_count is not None and len(read_count) == len(count_cols):
        weights = pd.Series(read_count, index=count_cols)
        weighted = merged[count_cols].multiply(weights, axis=1)
        return weighted.sum(axis=1) / weights.sum()
    else:
        # Fallback to simple average if weights are not provided
        return merged[count_cols].mean(axis=1)

def codon_weighted_average(merged: pd.DataFrame, codon_weights: Optional[List[float]] = None) -> pd.Series:
    """
    Calculate codon-weighted average across count columns.
    """
    count_cols = merged.columns[1:]
    if codon_weights is not None and len(codon_weights) == len(count_cols):
        weights = pd.Series(codon_weights, index=count_cols)
        weighted = merged[count_cols].multiply(weights, axis=1)
        return weighted.sum(axis=1) / weights.sum()
    else:
        # Fallback to simple average if weights are not provided
        return merged[count_cols].mean(axis=1)

def calculate_activity_scores(
    count_dfs: List[pd.DataFrame],
    method: str = 'rep-weighted',
    min_reads: int = 0,
    bins_required: int = 3,
    reps_required: int = 3,
    read_count: Optional[List[int]] = None,
    codon_weights: Optional[List[float]] = None
) -> pd.DataFrame:
    """
    Calculate activity scores for variants from count data.

    Parameters
    ----------
    count_dfs : list of pandas.DataFrame
        List of DataFrames, one per sample/replicate/bin, containing variant counts.
    method : str, default 'rep-weighted'
        Averaging method ('simple-avg', 'rep-weighted', 'codon-weighted').
    min_reads : int, default 0
        Minimum normalized reads per bin for a variant to be scored.
    bins_required : int, default 3
        Number of bins a variant must appear in to be scored.
    reps_required : int, default 3
        Number of replicates a variant must appear in to be scored.
    read_count : list of int, optional
        List of total reads per sample for normalization.
    codon_weights : list of float, optional
        List of weights for codon-weighted averaging.

    Returns
    -------
    scores_df : pandas.DataFrame
        DataFrame with activity scores for each variant.

    Examples
    --------
    >>> scores = calculate_activity_scores([df1, df2, df3], method='rep-weighted')
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting activity score calculation.")

    merged = merge_count_dfs(count_dfs)
    merged = normalize_counts(merged, read_count)
    merged = filter_variants(merged, min_reads, bins_required, reps_required)

    if method == 'simple-avg':
        merged['avgscore'] = simple_average(merged)
    elif method == 'rep-weighted':
        merged['avgscore_rep_weighted'] = rep_weighted_average(merged, read_count)
    elif method == 'codon-weighted':
        merged['avgscore_codon_weighted'] = codon_weighted_average(merged, codon_weights)
    else:
        logger.error(f"Unknown averaging method: {method}")
        raise ValueError(f"Unknown averaging method: {method}")

    logger.info("Activity score calculation complete.")
    return merged
