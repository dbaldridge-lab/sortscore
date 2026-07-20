"""
Data processing utilities for Sort-seq variant analysis.

This module provides functions for aggregating and transforming variant score data.

Examples
--------
>>> from sortscore.analysis.data_processing import aggregate_aa_data
>>> aa_data = aggregate_aa_data(scores_df, 'avgscore_rep_weighted')
"""
import pandas as pd


_AA_GROUP_COLUMNS = {'aa_seq_diff', 'annotate_aa'}
_AA_FIRST_COLUMNS = {
    'aa_seq',
    'batch',
}
_AA_SUM_COLUMNS = {'n_codons', 'n_measurements'}
_AA_DROPPED_COLUMNS = {
    'variant_seq',
    'codon_diff',
    'dna_seq_diff',
    'annotate_dna',
    'SD_codon',
    'SD_rep',
    'CV_rep',
    'CV_codon',
    'SEM',
    'CI_lower',
    'CI_upper',
}

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

def aggregate_synonymous_variants(scores_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate synonymous variants by averaging their scores.
    
    Groups variants by their AA sequence difference and variant annotation type 'synonymous',
    then averages individual replicate scores and recalculates avgscores.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with annotated scores including 'aa_seq_diff' and 'annotate_aa' columns.
        
    Returns
    -------
    aggregated_df : pd.DataFrame
        DataFrame with synonymous variants averaged into single rows.
        
    Examples
    --------
    >>> aggregated = aggregate_synonymous_variants(annotated_scores)
    """
    if 'aa_seq_diff' not in scores_df.columns or 'annotate_aa' not in scores_df.columns:
        raise ValueError("DataFrame must contain 'aa_seq_diff' and 'annotate_aa' columns")
    
    # Group by AA sequence difference and annotation type with proper aggregation for variance calculations
    # Use dropna=False to include NaN values (important for synonymous variants)
    agg_dict = {}
    
    for col in scores_df.columns:
        if col in _AA_GROUP_COLUMNS:
            continue
        if col.startswith('avgscore') or col.startswith('score.r') or (
            col.startswith('Rep') and col.endswith('.score')
        ):
            agg_dict[col] = 'mean'
        elif col.startswith('count.r') or col.startswith('norm.count.r') or (
            col.startswith('Rep') and col.endswith('.sum')
        ):
            agg_dict[col] = 'sum'
        elif col in _AA_SUM_COLUMNS:
            agg_dict[col] = 'sum'
        elif col in _AA_FIRST_COLUMNS:
            agg_dict[col] = 'first'
        elif col.startswith('prop.r') or col in _AA_DROPPED_COLUMNS:
            continue
        else:
            raise ValueError(f"No AA aggregation rule defined for column '{col}'.")
    
    aggregated_df = scores_df.groupby(['aa_seq_diff', 'annotate_aa'], dropna=False).agg(agg_dict).reset_index()
    
    return aggregated_df
