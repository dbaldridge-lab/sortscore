"""
Data processing utilities for Sort-seq variant analysis.

This module provides functions for aggregating and transforming variant score data.

Examples
--------
>>> from sortscore.analysis.data_processing import aggregate_aa_data
>>> aa_data = aggregate_aa_data(scores_df, 'avgscore_rep_weighted')
"""
import pandas as pd
from typing import Optional
from sortscore.sequence_parsing import translate_dna


def aggregate_aa_data(scores_df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """
    Aggregate DNA-level scores to amino acid level by averaging synonymous variants.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with DNA-level scores and aa_seq_diff column.
    score_col : str
        Name of the score column to aggregate.
    
    Returns
    -------
    aa_data : pd.DataFrame
        Aggregated AA-level data with averaged scores.
        
    Examples
    --------
    >>> aa_data = aggregate_aa_data(dna_scores, 'avgscore_rep_weighted')
    """
    if 'aa_seq_diff' not in scores_df.columns:
        raise ValueError("scores_df must contain 'aa_seq_diff' column for AA aggregation")
    
    # Group by aa_seq_diff and aggregate scores
    agg_cols = {score_col: 'mean'}
    if 'annotate_aa' in scores_df.columns:
        agg_cols['annotate_aa'] = 'first'
    
    aa_data = scores_df.groupby('aa_seq_diff').agg(agg_cols).reset_index()
    
    return aa_data
