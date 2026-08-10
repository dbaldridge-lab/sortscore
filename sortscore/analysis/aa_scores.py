"""
Amino acid scores processing and export for Sort-seq analysis.

This module provides functions for processing amino acid scores from DNA variant data,
including aggregation of synonymous codons and file export.
"""
import os
import logging
import pandas as pd
from typing import List


def _get_score_column_from_avg_method(avg_method: str) -> str:
    """Return the canonical output column for a configured averaging method."""
    if avg_method not in {'simple-avg', 'rep-weighted'}:
        raise ValueError(f"Unknown averaging method: {avg_method}")
    return 'score'


def _get_count_columns(scores_df: pd.DataFrame) -> List[str]:
    """Return raw count columns that should be preserved in AA score outputs."""
    return [
        col for col in scores_df.columns
        if col.startswith('count.r') and 'b' in col
    ]


def process_and_save_aa_scores(scores_df: pd.DataFrame, experiment, scores_dir: str, analysis_logger) -> None:
    """
    Process and save amino acid scores from variant data.
    
    This function handles the complete AA scores workflow including:
    - Filtering out NaN values
    - Checking if codon aggregation is needed
    - Saving to CSV file
    - Logging output
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing variant scores and annotations
    experiment : ExperimentConfig
        Experiment configuration containing metadata
    scores_dir : str
        Directory to save scores file
    analysis_logger : AnalysisLogger
        Logger instance for recording outputs
        
    Examples
    --------
    >>> process_and_save_aa_scores(scores_df, experiment, 'output/scores', 'suffix', logger)
    """
    if 'aa_seq_diff' not in scores_df.columns:
        return

    score_col = _get_score_column_from_avg_method(experiment.avg_method)
    
    aa_scores = build_aa_scores_table(scores_df, score_col)
    
    # Save to file
    aa_scores_file = os.path.join(scores_dir, f"{experiment.experiment_name}_aa_scores.csv")
    aa_scores.to_csv(aa_scores_file, index=False)
    logging.info(f"Saved AA scores to {aa_scores_file} ({len(aa_scores)} unique AA variants)")
    
    # Log file output  
    analysis_logger.log_output_file(
        'aa_scores',
        f"{experiment.experiment_name}_aa_scores.csv",
        aa_scores_file,
        variant_count=len(aa_scores)
    )


def build_aa_scores_table(scores_df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """
    Build the AA-level score table from a scored variant table.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing variant scores and annotations.
    score_col : str
        Score column used to determine which rows are retained and how codon
        variance is calculated.

    Returns
    -------
    pd.DataFrame
        AA-level score table with synonymous codons aggregated when needed.
    """
    if 'aa_seq_diff' not in scores_df.columns:
        raise ValueError("scores_df must contain 'aa_seq_diff' column for AA score output")

    scores_df_drop_nan = scores_df.dropna(subset=[score_col])
    rep_score_columns = [
        col for col in scores_df_drop_nan.columns
        if col.startswith('Rep') and col.endswith('.score')
    ]
    return _check_codon_num(scores_df_drop_nan, score_col, rep_score_columns)


def _check_codon_num(scores_df_drop_nan: pd.DataFrame, score_col: str, 
                                  rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Check if codon aggregation is needed and process AA scores accordingly.
    
    This function checks if there are multiple codons per AA variant and processes
    the data using either DNA->AA aggregation or AA-only score export.
    
    Parameters
    ----------
    scores_df_drop_nan : pd.DataFrame
        DataFrame with NaN values already filtered out
    score_col : str
        Name of the score column to use
    rep_score_columns : List[str]
        List of replicate score column names
        
    Returns
    -------
    pd.DataFrame
        Processed AA scores
    """
    # Check if there are multiple codons per AA variant (DNA->AA case)
    aa_variant_counts = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa']).size()
    needs_aggregation = (aa_variant_counts > 1).any()
    
    if needs_aggregation:
        return _process_dna_to_aa_aggregation(scores_df_drop_nan, score_col, rep_score_columns)
    else:
        return _process_aa_only_scores(scores_df_drop_nan, rep_score_columns)


def _process_dna_to_aa_aggregation(scores_df_drop_nan: pd.DataFrame, score_col: str, 
                                  rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Process DNA->AA aggregation case.
    
    Parameters
    ----------
    scores_df_drop_nan : pd.DataFrame
        DataFrame with variant scores (NaN filtered)
    score_col : str
        Name of the score column to use
    rep_score_columns : List[str]
        List of replicate score column names
        
    Returns
    -------
    pd.DataFrame
        Aggregated AA scores
    """
    # DNA->AA aggregation case: aggregate synonymous variants
    count_columns = _get_count_columns(scores_df_drop_nan)
    columns_to_average = [score_col] + rep_score_columns
    
    # Calculate standard deviation and count of codon-level scores before AA aggregation
    aa_scores_std = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[score_col].agg(['std', 'count']).reset_index()
    aa_scores_std.columns = ['aa_seq_diff', 'annotate_aa', 'SD_codon', 'n_codons']
    
    # Calculate mean scores for aggregation
    aa_scores = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[columns_to_average].mean().reset_index()

    if count_columns:
        aa_count_sums = (
            scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[count_columns]
            .sum()
            .reset_index()
        )
        aa_scores = aa_scores.merge(aa_count_sums, on=['aa_seq_diff', 'annotate_aa'], how='left')
    
    # Preserve codon-level score spread/counts for aggregated AA rows.
    aa_scores = aa_scores.merge(aa_scores_std, on=['aa_seq_diff', 'annotate_aa'], how='left')

    return aa_scores


def _process_aa_only_scores(scores_df_drop_nan: pd.DataFrame, rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Process AA-only case.
    
    Parameters
    ----------
    scores_df_drop_nan : pd.DataFrame
        DataFrame with variant scores (NaN filtered)
    rep_score_columns : List[str]
        List of replicate score column names
        
    Returns
    -------
    pd.DataFrame
        AA scores
    """
    # AA-only case: no aggregation needed, just copy the data
    count_columns = _get_count_columns(scores_df_drop_nan)
    columns_to_include = ['aa_seq_diff', 'annotate_aa', 'score'] + rep_score_columns + count_columns
    
    return scores_df_drop_nan[columns_to_include].copy()
