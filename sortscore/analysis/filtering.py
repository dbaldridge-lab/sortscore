"""
Filtering utilities for Sort-seq variant analysis.

This module provides functions to filter variants by minimum reads, bins, and replicates required.

Examples
--------
>>> from sortscore.analysis.filtering import filter_variants
"""
import pandas as pd
from typing import Optional

def filter_variants(
    merged: pd.DataFrame,
    min_reads: int,
    bins_required: int,
    reps_required: int
) -> pd.DataFrame:
    """
    Filter variants by minimum reads, bins, and replicates required.

    Parameters
    ----------
    merged : pandas.DataFrame
        DataFrame with variant counts (wide format).
    min_reads : int
        Minimum normalized reads per bin for a variant to be scored.
    bins_required : int
        Number of bins a variant must appear in to be scored.
    reps_required : int
        Number of replicates a variant must appear in to be scored.

    Returns
    -------
    filtered : pandas.DataFrame
        Filtered DataFrame.

    Examples
    --------
    >>> filtered = filter_variants(df, 0, 3, 3)
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
