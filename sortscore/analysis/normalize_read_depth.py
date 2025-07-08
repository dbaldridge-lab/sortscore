"""
Normalization utilities for Sort-seq variant analysis.

This module provides functions to normalize variant counts per million reads and apply minimum read thresholds, as in the original notebook.

Examples
--------
>>> from sortscore.analysis.normalization import normalize_counts_per_million, apply_minread_threshold
>>> df_norm = normalize_counts_per_million(df, read_count)
>>> df_norm = apply_minread_threshold(df_norm, minread_threshold)
"""
import pandas as pd
import numpy as np
from typing import List, Dict

def counts_per_million(df: pd.DataFrame, count_cols: List[str], read_counts: List[float]) -> pd.DataFrame:
    """
    Normalize variant counts to reads per million for each count column.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing raw counts.
    count_cols : list of str
        List of column names with raw counts.
    read_counts : list of float
        Total read counts for each column (same order as count_cols).

    Returns
    -------
    df_norm : pandas.DataFrame
        DataFrame with normalized count columns (named 'norm.' + col).
    """
    df_norm = df.copy()
    for col, rc in zip(count_cols, read_counts):
        norm_col = f'norm.{col}'
        df_norm[norm_col] = df_norm[col] / rc * 1e6
    return df_norm

def apply_minread_threshold(df: pd.DataFrame, norm_cols: List[str], minread_threshold: float) -> pd.DataFrame:
    """
    Apply a minimum read threshold to normalized count columns, setting values below threshold to NaN.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with normalized count columns.
    norm_cols : list of str
        List of normalized count column names.
    minread_threshold : float
        Minimum reads per million required to keep a value.

    Returns
    -------
    df_out : pandas.DataFrame
        DataFrame with values below threshold set to NaN in normalized columns.
    """
    df_out = df.copy()
    for col in norm_cols:
        df_out[col] = df_out[col].where(df_out[col] >= minread_threshold, np.nan)
    return df_out

"""
Proportion calculation utilities for Sort-seq variant analysis.

This module provides a function to calculate the proportion of reads per bin for each variant, as in the notebook.

Examples
--------
>>> from sortscore.analysis.normalization import calculate_bin_proportions
>>> df = calculate_bin_proportions(df, norm_cols)
"""
import pandas as pd
from typing import List

def calculate_bin_proportions(df: pd.DataFrame, norm_cols: List[str], prop_prefix: str = 'prop') -> pd.DataFrame:
    """
    Calculate the proportion of normalized reads in each bin relative to the total normalized reads for that replicate (row-wise).

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with normalized count columns (e.g., 'norm.counts.2', ...).
    norm_cols : list of str
        List of normalized count column names (one per bin, for a single replicate).
    prop_prefix : str, default 'prop'
        Prefix for the output proportion columns.

    Returns
    -------
    df_out : pandas.DataFrame
        DataFrame with new columns for the proportion of reads in each bin (for that replicate).
    """
    df_out = df.copy()
    total = df_out[norm_cols].sum(axis=1)
    for col in norm_cols:
        prop_col = f'{prop_prefix}.{col}'
        df_out[prop_col] = df_out[col] / total
    return df_out
