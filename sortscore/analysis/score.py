"""
Score calculation functions for variant activity analysis.

This module provides functions to calculate activity scores from count data, supporting different averaging methods and filtering criteria.

Examples
--------
>>> from sortscore.analysis.score import calculate_activity_scores
"""
import logging
import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional
from sortscore.analysis.filtering import filter_variants

def filter_by_cv(df: pd.DataFrame, score_cols: List[str], max_cv: float) -> pd.Series:
    """
    Filter variants by coefficient of variation across replicates.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing replicate scores
    score_cols : list of str
        Column names containing replicate scores
    max_cv : float
        Maximum CV threshold as decimal (0.5 = 50% CV)
        Variants with CV > max_cv will be filtered out (quality filter)
        
    Returns
    -------
    pd.Series
        Boolean mask for variants that pass CV filter
    """
    cv = df[score_cols].std(axis=1) / df[score_cols].mean(axis=1)
    return (cv <= max_cv) | cv.isna()

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


def calculate_activity_scores(
    count_dfs: List[pd.DataFrame],
    method: str = 'rep-weighted',
    min_reads: int = 0,
    bins_required: int = 1,
    reps_required: int = 3,
    read_count: Optional[List[int]] = None
) -> pd.DataFrame:
    """
    Calculate activity scores for variants from count data.

    Parameters
    ----------
    count_dfs : list of pandas.DataFrame
        List of DataFrames, one per sample/replicate/bin, containing variant counts.
    method : str, default 'rep-weighted'
        Averaging method ('simple-avg', 'rep-weighted').
    min_reads : int, default 0
        Minimum normalized reads per bin for a variant to be scored.
    bins_required : int, default 1
        Number of bins a variant must appear in to be scored.
    reps_required : int, default 3
        Number of replicates a variant must appear in to be scored.
    read_count : list of int, optional
        List of total reads per sample for normalization.

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
    else:
        logger.error(f"Unknown averaging method: {method}")
        raise ValueError(f"Unknown averaging method: {method}")

    logger.info("Activity score calculation complete.")
    return merged

def calculate_full_activity_scores(
    counts: Dict[int, Dict[int, pd.DataFrame]],
    median_gfp: Dict[int, Dict[int, float]],
    merged_df: pd.DataFrame,
    min_bins: int = 3,
    min_reps: int = 3,
    minread_threshold: float = 0.0,
    avg_method: str = 'rep-weighted',
    groupby_cols: Optional[list] = None,
    total_reads: Optional[Dict[int, Dict[int, int]]] = None,
    cell_prop: Optional[Dict[int, Dict[int, float]]] = None,
    max_cv: Optional[float] = None,
) -> pd.DataFrame:
    """
    Calculate activity scores for all variants using full Sort-seq logic (per-bin/rep normalization, bin proportions, replicate/rep-weighted averaging).

    Parameters
    ----------
    counts : dict
        Nested dict of DataFrames: counts[rep][bin] = DataFrame with columns ['variant_seq', 'count'] (or similar).
    median_gfp : dict
        Nested dict of median GFP values: median_gfp[rep][bin] = float.
    min_bins : int, default 3
        Minimum number of bins a variant must appear in per replicate to be scored.
    min_reps : int, default 3
        Minimum number of replicates a variant must appear in to be scored.
    minread_threshold : float, default 0.0
        Minimum normalized reads per million required to keep a value.
    avg_method : str, default 'rep-weighted'
        Averaging method: 'simple-avg' or 'rep-weighted'.
    groupby_cols : list, optional
        Columns to group by for codon/synonymous averaging (e.g., ['annotate_aa', 'aa_seq_diff']).
    total_reads : dict, optional
        Nested dict of total sequencing reads for normalization: total_reads[rep][bin] = int.
        If not provided, uses sum of variant counts in each bin.
    cell_prop : dict, optional
        Nested dict of cell proportions per bin for normalization: cell_prop[rep][bin] = float.
        Values should be proportions (0.0-1.0) representing the %Gate values from cell sorter output.
        Proportions do NOT need to sum to 1.0 across bins.
        If provided, normalizes by both sequencing depth and cell proportions.

    Returns
    -------
    scores_df : pd.DataFrame
        DataFrame with all replicate/bin scores, averages, and annotations.

    Examples
    --------
    >>> scores = calculate_full_activity_scores(counts, median_gfp)
    """
    # 1. Use the provided merged DataFrame
    df = merged_df.copy()

    # 3. Normalize counts per million for each rep/bin (accounting for sequencing depth and cell proportions)
    for rep in counts:
        for bin in counts[rep]:
            col = f'count.r{rep}b{bin}'
            # Use external total reads (prior to filtering) if provided, otherwise sum of variant counts
            if total_reads is not None and rep in total_reads and bin in total_reads[rep]:
                total_reads_for_norm = total_reads[rep][bin]
            else:
                total_reads_for_norm = counts[rep][bin]['count'].sum()
            
            norm_col = f'norm.{col}'
            # Base normalization: reads per million
            df[norm_col] = df[col] / total_reads_for_norm * 1e6
            
            # Additional normalization by cell proportions if available
            if cell_prop is not None and rep in cell_prop and bin in cell_prop[rep]:
                cell_proportion = cell_prop[rep][bin]
                if cell_proportion > 0:  # Avoid division by zero
                    df[norm_col] = df[norm_col] / cell_proportion
            
            # Apply minread threshold
            if minread_threshold > 0:
                df[norm_col] = df[norm_col].where(df[norm_col] >= minread_threshold, np.nan)

    # 4. For each replicate, sum normalized counts across bins, require min_bins
    rep_sums = {}
    for rep in counts:
        norm_cols = [f'norm.count.r{rep}b{bin}' for bin in counts[rep]]
        sum_col = f'Rep{rep}.sum'
        df[sum_col] = df[norm_cols].sum(axis=1, min_count=min_bins)
        rep_sums[rep] = sum_col

    # 5. Calculate bin proportions for each rep/bin
    for rep in counts:
        norm_cols = [f'norm.count.r{rep}b{bin}' for bin in counts[rep]]
        for i, bin in enumerate(counts[rep]):
            prop_col = f'prop.r{rep}b{bin}'
            df[prop_col] = df[f'norm.count.r{rep}b{bin}'] / df[f'Rep{rep}.sum']

    # 6. Calculate bin activity scores (proportion * median GFP)
    for rep in counts:
        for bin in counts[rep]:
            score_col = f'score.r{rep}b{bin}'
            prop_col = f'prop.r{rep}b{bin}'
            gfp = median_gfp[rep][bin]
            df[score_col] = df[prop_col] * gfp

    # 7. Replicate-level activity score: sum bin scores, require min_bins
    for rep in counts:
        score_cols = [f'score.r{rep}b{bin}' for bin in counts[rep]]
        rep_score_col = f'Rep{rep}.score'
        df[rep_score_col] = df[score_cols].sum(axis=1, min_count=min_bins)

    # 8. Average scores (simple, rep-weighted)
    rep_score_cols = [f'Rep{rep}.score' for rep in counts]
    # Only keep rows with at least min_reps non-NaN replicate scores
    valid = df[rep_score_cols].notna().sum(axis=1) >= min_reps
    
    # Apply CV filtering if specified
    if max_cv is not None:
        cv_valid = filter_by_cv(df, rep_score_cols, max_cv)
        valid = valid & cv_valid
    
    df_valid = df[valid].copy()
    # Simple average
    df_valid['avgscore'] = df_valid[rep_score_cols].mean(axis=1)
    # Replicate-weighted average
    total_weight = sum(df_valid[f'Rep{rep}.sum'].fillna(0) for rep in counts)
    rep_weighted_score = sum(df_valid[f'Rep{rep}.score'].fillna(0) * df_valid[f'Rep{rep}.sum'].fillna(0) for rep in counts)
    df_valid['avgscore_rep_weighted'] = rep_weighted_score / total_weight

    # 10. Optionally, group by synonymous/codon group for AA-level scores
    if groupby_cols is not None:
        group = df_valid.groupby(groupby_cols)
        df_aa = group[[
            'avgscore',
            'avgscore_rep_weighted',
            *rep_score_cols
        ]].mean().reset_index()
        return df_valid, df_aa
    return df_valid
