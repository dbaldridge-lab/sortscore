"""
Statistical analysis functions for sortscore data processing.

This module provides functions for calculating replicate statistics, confidence intervals,
and other statistical measures commonly used in Sort-seq variant analysis.
"""
import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from typing import List


def calculate_replicate_statistics(df: pd.DataFrame, rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Calculate statistical measures for replicate scores including SD, CV, and 95% CI.
    
    This function computes standard deviation, coefficient of variation, and 95% confidence
    intervals for replicate measurements using t-distribution with appropriate degrees of freedom.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing replicate score columns
    rep_score_columns : List[str]
        List of column names containing replicate scores (e.g., ['Rep1.score', 'Rep2.score'])
        
    Returns
    -------
    pd.DataFrame
        DataFrame with added statistical columns:
        - SD_rep: Standard deviation of replicates (rounded to integers)
        - CV_rep: Coefficient of variation (3 decimal places)
        - n_measurements: Number of non-null measurements per variant
        - CI_lower: Lower bound of 95% confidence interval
        - CI_upper: Upper bound of 95% confidence interval
        
    Examples
    --------
    >>> rep_cols = ['Rep1.score', 'Rep2.score', 'Rep3.score']
    >>> df_with_stats = calculate_replicate_statistics(df, rep_cols)
    """
    df_result = df.copy()
    
    if len(rep_score_columns) < 2:
        return df_result
    
    # Calculate basic statistics
    rep_mean = df_result[rep_score_columns].mean(axis=1)
    rep_std = df_result[rep_score_columns].std(axis=1, ddof=1)
    
    # Calculate n_measurements dynamically based on non-empty replicate values
    # Each row represents one variant, so n_measurements = number of non-null replicates
    n_measurements = df_result[rep_score_columns].notna().sum(axis=1)
    
    # Calculate 95% CI using t-distribution with actual degrees of freedom
    df_actual = n_measurements - 1
    t_critical = scipy_stats.t.ppf(0.975, df_actual)
    margin_of_error = t_critical * (rep_std / np.sqrt(n_measurements))
    
    # Add stats to DataFrame
    df_result['SD_rep'] = rep_std.round().astype('Int64')
    df_result['CV_rep'] = (rep_std / rep_mean).round(3)
    df_result['n_measurements'] = n_measurements.astype('Int64')
    df_result['CI_lower'] = (rep_mean - margin_of_error).round().astype('Int64')
    df_result['CI_upper'] = (rep_mean + margin_of_error).round().astype('Int64')
    
    return df_result


def round_score_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Round score columns to integers for cleaner output.
    
    Identifies columns containing 'score' in the name and rounds floating-point
    values to integers using pandas nullable integer type.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing score columns to round
        
    Returns
    -------
    pd.DataFrame
        DataFrame with score columns rounded to integers
        
    Examples
    --------
    >>> df_rounded = round_score_columns(df)
    """
    df_result = df.copy()
    
    # Find and round score columns
    score_columns = [col for col in df_result.columns if 'score' in col.lower()]
    for col in score_columns:
        if df_result[col].dtype in ['float64', 'float32']:
            df_result[col] = df_result[col].round().astype('Int64')
    
    return df_result


def get_replicate_score_columns(df: pd.DataFrame) -> List[str]:
    """
    Identify replicate score columns in a DataFrame.
    
    Looks for columns matching the pattern 'Rep*.score' (e.g., Rep1.score, Rep2.score).
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to search for replicate columns
        
    Returns
    -------
    List[str]
        List of replicate score column names
        
    Examples
    --------
    >>> rep_cols = get_replicate_score_columns(df)
    >>> print(rep_cols)
    ['Rep1.score', 'Rep2.score', 'Rep3.score']
    """
    return [col for col in df.columns if col.startswith('Rep') and col.endswith('.score')]


def calculate_codon_and_replicate_variance(aa_scores: pd.DataFrame, rep_score_columns: List[str]) -> pd.DataFrame:
    """
    Calculate statistics with codon and replicate variance decomposition.
    
    This function computes advanced statistics for amino acid aggregated data that includes
    both codon-level variance (from multiple codons per amino acid) and replicate variance.
    It calculates total variance as the sum of codon and replicate variances, then computes
    standard error and confidence intervals considering the total sample size.
    
    Parameters
    ----------
    aa_scores : pd.DataFrame
        DataFrame containing aggregated amino acid scores with 'SD_codon' and 'n_codons' columns
    rep_score_columns : List[str]
        List of replicate score column names
        
    Returns
    -------
    pd.DataFrame
        DataFrame with added statistical columns:
        - SD_rep: Standard deviation of replicates
        - CV_rep: Coefficient of variation for replicates
        - CV_codon: Coefficient of variation for codons
        - n_measurements: Total number of measurements (codons × replicates)
        - SEM: Standard error of the mean
        - CI_lower: Lower bound of 95% confidence interval
        - CI_upper: Upper bound of 95% confidence interval
        - SD_codon: Rounded to integers
        
    Examples
    --------
    >>> rep_cols = ['Rep1.score', 'Rep2.score', 'Rep3.score']
    >>> aa_scores_with_stats = calculate_codon_and_replicate_variance(aa_scores, rep_cols)
    """
    df_result = aa_scores.copy()
    
    if len(rep_score_columns) < 2:
        return df_result
    
    # Calculate replicate statistics
    aa_rep_mean = df_result[rep_score_columns].mean(axis=1)
    aa_rep_std = df_result[rep_score_columns].std(axis=1, ddof=1)
    
    # Calculate variances
    var_codon = df_result['SD_codon'] ** 2
    var_rep = aa_rep_std ** 2
    
    # Calculate total variance (sum of codon and replicate variances)
    var_total = var_codon.fillna(0) + var_rep
    
    # Calculate n_measurements dynamically
    n_measurements = df_result['n_codons'] * df_result[rep_score_columns].notna().sum(axis=1)
    
    # Calculate SEM considering total sample size
    sem = np.sqrt(var_total / n_measurements)
    
    # Calculate 95% CI using t-distribution
    df_total = n_measurements - 1
    t_critical = scipy_stats.t.ppf(0.975, df_total)
    aa_margin_of_error = t_critical * sem
    
    # Add statistics to DataFrame
    df_result['SD_rep'] = aa_rep_std.round().astype('Int64')
    df_result['CV_rep'] = (aa_rep_std / aa_rep_mean).round(3)
    df_result['CV_codon'] = (df_result['SD_codon'] / aa_rep_mean).round(3)
    df_result['n_measurements'] = n_measurements.astype('Int64')
    df_result['SEM'] = sem.round().astype('Int64')
    df_result['CI_lower'] = (aa_rep_mean - aa_margin_of_error).round().astype('Int64')
    df_result['CI_upper'] = (aa_rep_mean + aa_margin_of_error).round().astype('Int64')
    
    # Round SD_codon column to integers
    if 'SD_codon' in df_result.columns:
        df_result['SD_codon'] = df_result['SD_codon'].round().astype('Int64')
    
    return df_result