"""
Amino acid scores processing and export for Sort-seq analysis.

This module provides functions for processing amino acid scores from DNA variant data,
including aggregation of synonymous codons, statistical analysis, and file export.
"""
import os
import logging
import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from typing import Tuple, List
from sortscore.analysis.statistics import calculate_codon_and_replicate_variance


def process_and_save_aa_scores(scores_df: pd.DataFrame, experiment, scores_dir: str, 
                              output_suffix: str, analysis_logger) -> None:
    """
    Process and save amino acid scores from variant data.
    
    This function handles the complete AA scores workflow including:
    - Filtering out NaN values
    - Checking if codon aggregation is needed
    - Calculating appropriate statistics (with/without codon variance)
    - Rounding score columns
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
    output_suffix : str
        Suffix for output filename
    analysis_logger : AnalysisLogger
        Logger instance for recording outputs
        
    Examples
    --------
    >>> process_and_save_aa_scores(scores_df, experiment, 'output/scores', 'suffix', logger)
    """
    if 'aa_seq_diff' not in scores_df.columns:
        return
        
    # Determine score column
    if experiment.avg_method == 'simple-avg':
        score_col = 'avgscore'
    else:
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
    
    # Filter out rows with NaN values first
    scores_df_drop_nan = scores_df.dropna(subset=[score_col])
    
    # Find replicate score columns
    rep_score_columns = [col for col in scores_df_drop_nan.columns 
                        if col.startswith('Rep') and col.endswith('.score')]
    
    # Check aggregation needs and process scores
    aa_scores = _check_codon_num(scores_df_drop_nan, score_col, rep_score_columns)
    
    # Round score columns to integers
    aa_scores = _round_score_columns(aa_scores)
    
    # Save to file
    aa_scores_file = os.path.join(scores_dir, f"{experiment.experiment_name}_aa_scores_{output_suffix}.csv")
    aa_scores.to_csv(aa_scores_file, index=False)
    logging.info(f"Saved AA scores to {aa_scores_file} ({len(aa_scores)} unique AA variants)")
    
    # Log file output  
    analysis_logger.log_output_file(
        'aa_scores',
        f"{experiment.experiment_name}_aa_scores_{output_suffix}.csv", 
        aa_scores_file,
        variant_count=len(aa_scores)
    )


def _check_codon_num(scores_df_drop_nan: pd.DataFrame, score_col: str, 
                                  rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Check if codon aggregation is needed and process AA scores accordingly.
    
    This function checks if there are multiple codons per AA variant and processes
    the data using either DNA->AA aggregation (with codon variance) or AA-only
    statistics (replicate variance only).
    
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
        Processed AA scores with appropriate statistics
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
    Process DNA->AA aggregation case with codon variance decomposition.
    
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
        Aggregated AA scores with codon and replicate statistics
    """
    # DNA->AA aggregation case: aggregate synonymous variants
    columns_to_average = ['avgscore', 'avgscore_rep_weighted'] + rep_score_columns
    
    # Calculate standard deviation and count of codon-level scores before AA aggregation
    aa_scores_std = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[score_col].agg(['std', 'count']).reset_index()
    aa_scores_std.columns = ['aa_seq_diff', 'annotate_aa', 'SD_codon', 'n_codons']
    
    # Calculate mean scores for aggregation
    aa_scores = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[columns_to_average].mean().reset_index()
    
    # Merge the standard deviation and count of codon scores
    aa_scores = aa_scores.merge(aa_scores_std, on=['aa_seq_diff', 'annotate_aa'], how='left')
    
    # Calculate statistics with codon and replicate variance decomposition
    aa_scores = calculate_codon_and_replicate_variance(aa_scores, rep_score_columns)
    
    return aa_scores


def _process_aa_only_scores(scores_df_drop_nan: pd.DataFrame, rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Process AA-only case with simple replicate statistics.
    
    Parameters
    ----------
    scores_df_drop_nan : pd.DataFrame
        DataFrame with variant scores (NaN filtered)
    rep_score_columns : List[str]
        List of replicate score column names
        
    Returns
    -------
    pd.DataFrame
        AA scores with replicate statistics only
    """
    # AA-only case: no aggregation needed, just copy the data
    columns_to_include = ['aa_seq_diff', 'annotate_aa', 'avgscore', 'avgscore_rep_weighted'] + rep_score_columns
    
    aa_scores = scores_df_drop_nan[columns_to_include].copy()
    
    # Calculate simple replicate statistics (no codon variance)
    if len(rep_score_columns) >= 2:
        aa_rep_mean = aa_scores[rep_score_columns].mean(axis=1)
        aa_rep_std = aa_scores[rep_score_columns].std(axis=1, ddof=1)
        
        # Calculate n_measurements (just number of non-empty replicates)
        n_measurements = aa_scores[rep_score_columns].notna().sum(axis=1)
        
        # Calculate SEM using only replicate variance
        sem = aa_rep_std / np.sqrt(n_measurements)
        
        # Calculate 95% CI using t-distribution
        df_actual = n_measurements - 1
        t_critical = scipy_stats.t.ppf(0.975, df_actual)
        aa_margin_of_error = t_critical * sem
        
        aa_scores['SD_rep'] = aa_rep_std.round().astype('Int64')
        aa_scores['CV_rep'] = (aa_rep_std / aa_rep_mean).round(3)
        aa_scores['n_measurements'] = n_measurements.astype('Int64')
        aa_scores['SEM'] = sem.round().astype('Int64')
        aa_scores['CI_lower'] = (aa_rep_mean - aa_margin_of_error).round().astype('Int64')
        aa_scores['CI_upper'] = (aa_rep_mean + aa_margin_of_error).round().astype('Int64')
    
    return aa_scores


def _round_score_columns(aa_scores: pd.DataFrame) -> pd.DataFrame:
    """
    Round score columns to integers for cleaner output.
    
    Parameters
    ----------
    aa_scores : pd.DataFrame
        DataFrame containing score columns
        
    Returns
    -------
    pd.DataFrame
        DataFrame with score columns rounded to integers
    """
    # Round score columns to integers
    score_columns = [col for col in aa_scores.columns if 'score' in col.lower()]
    for col in score_columns:
        if aa_scores[col].dtype in ['float64', 'float32']:
            aa_scores[col] = aa_scores[col].round().astype('Int64')
    
    return aa_scores