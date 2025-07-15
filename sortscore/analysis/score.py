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

def calculate_full_activity_scores(
    counts: Dict[int, Dict[int, pd.DataFrame]],
    median_gfp: Dict[int, Dict[int, float]],
    min_bins: int = 3,
    min_reps: int = 3,
    minread_threshold: float = 0.0,
    avg_method: str = 'rep-weighted',
    groupby_cols: Optional[list] = None,
    total_reads: Optional[Dict[int, Dict[int, int]]] = None,
) -> pd.DataFrame:
    """
    Calculate activity scores for all variants using full Sort-seq logic (per-bin/rep normalization, bin proportions, replicate/codon/rep-weighted averaging).

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
        Averaging method: 'simple-avg', 'rep-weighted', or 'codon-weighted'.
    groupby_cols : list, optional
        Columns to group by for codon/synonymous averaging (e.g., ['annotate_aa', 'aa_seq_diff']).
    total_reads : dict, optional
        Nested dict of total sequencing reads for normalization: total_reads[rep][bin] = int.
        If not provided, uses sum of variant counts in each bin.

    Returns
    -------
    scores_df : pd.DataFrame
        DataFrame with all replicate/bin scores, averages, and annotations.

    Examples
    --------
    >>> scores = calculate_full_activity_scores(counts, median_gfp)
    """
    # 1. Merge all counts into a single DataFrame (wide format)
    all_variants = set()
    for rep in counts:
        for bin_ in counts[rep]:
            all_variants.update(counts[rep][bin_]['variant_seq'])
    all_variants = sorted(all_variants)
    df = pd.DataFrame({'variant_seq': all_variants})

    # 2. Add count columns for each rep/bin
    for rep in counts:
        for bin_ in counts[rep]:
            col = f'count.r{rep}b{bin_}'
            d = counts[rep][bin_].set_index('variant_seq')['count']
            df[col] = df['variant_seq'].map(d).fillna(0)

    # 3. Normalize counts per million for each rep/bin
    for rep in counts:
        for bin_ in counts[rep]:
            col = f'count.r{rep}b{bin_}'
            # Use external total reads (prior to filtering) if provided, otherwise sum of variant counts
            if total_reads is not None and rep in total_reads and bin_ in total_reads[rep]:
                total_reads_for_norm = total_reads[rep][bin_]
            else:
                total_reads_for_norm = counts[rep][bin_]['count'].sum()
            norm_col = f'norm.{col}'
            df[norm_col] = df[col] / total_reads_for_norm * 1e6
            # Apply minread threshold
            if minread_threshold > 0:
                df[norm_col] = df[norm_col].where(df[norm_col] >= minread_threshold, np.nan)

    # 4. For each replicate, sum normalized counts across bins, require min_bins
    rep_sums = {}
    for rep in counts:
        norm_cols = [f'norm.count.r{rep}b{bin_}' for bin_ in counts[rep]]
        sum_col = f'Rep{rep}.sum'
        df[sum_col] = df[norm_cols].sum(axis=1, min_count=min_bins)
        rep_sums[rep] = sum_col

    # 5. Calculate bin proportions for each rep/bin
    for rep in counts:
        norm_cols = [f'norm.count.r{rep}b{bin_}' for bin_ in counts[rep]]
        for i, bin_ in enumerate(counts[rep]):
            prop_col = f'prop.r{rep}b{bin_}'
            df[prop_col] = df[f'norm.count.r{rep}b{bin_}'] / df[f'Rep{rep}.sum']

    # 6. Calculate bin activity scores (proportion * median GFP)
    for rep in counts:
        for bin_ in counts[rep]:
            score_col = f'score.r{rep}b{bin_}'
            prop_col = f'prop.r{rep}b{bin_}'
            gfp = median_gfp[rep][bin_]
            df[score_col] = df[prop_col] * gfp

    # 7. Replicate-level activity score: sum bin scores, require min_bins
    for rep in counts:
        score_cols = [f'score.r{rep}b{bin_}' for bin_ in counts[rep]]
        rep_score_col = f'Rep{rep}.score'
        df[rep_score_col] = df[score_cols].sum(axis=1, min_count=min_bins)

    # 8. Replicate/codon weights for codon-weighted averaging
    # (For now, use sum of RepX.sum as weights; can be extended for synonymous/codon groups)
    df['total_syn'] = sum(df[f'Rep{rep}.sum'].fillna(0) for rep in counts)
    for rep in counts:
        df[f'Rep{rep}.cw'] = df[f'Rep{rep}.sum'] / df['total_syn']
        df[f'Rep{rep}.score.cw'] = df[f'Rep{rep}.score'] * df[f'Rep{rep}.cw']

    # 9. Average scores (simple, rep-weighted, codon-weighted)
    rep_score_cols = [f'Rep{rep}.score' for rep in counts]
    rep_score_cw_cols = [f'Rep{rep}.score.cw' for rep in counts]
    # Only keep rows with at least min_reps non-NaN replicate scores
    valid = df[rep_score_cols].notna().sum(axis=1) >= min_reps
    df_valid = df[valid].copy()
    # Simple average
    df_valid['avgscore'] = df_valid[rep_score_cols].mean(axis=1)
    # Replicate-weighted average
    total_weight = sum(df_valid[f'Rep{rep}.sum'].fillna(0) for rep in counts)
    rep_weighted_score = sum(df_valid[f'Rep{rep}.score'].fillna(0) * df_valid[f'Rep{rep}.sum'].fillna(0) for rep in counts)
    df_valid['avgscore_rep_weighted'] = rep_weighted_score / total_weight
    # Codon-weighted average (sum of RepX.score.cw)
    df_valid['avgscore_codon_weighted'] = df_valid[rep_score_cw_cols].sum(axis=1)

    # 10. Optionally, group by synonymous/codon group for AA-level scores
    if groupby_cols is not None:
        group = df_valid.groupby(groupby_cols)
        df_aa = group[[
            'avgscore',
            'avgscore_rep_weighted',
            'avgscore_codon_weighted',
            *rep_score_cols
        ]].mean().reset_index()
        return df_valid, df_aa
    return df_valid
