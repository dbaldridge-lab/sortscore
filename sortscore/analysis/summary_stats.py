"""
Summary stats calculation and export for Sort-seq analysis.

This module provides functions for calculating summary stats
from Sort-seq variant analysis results and saving them to JSON format.
"""
import json
import os
import logging
import pandas as pd
from typing import Dict, Any, Optional
from sortscore.analysis.variant_aggregation import aggregate_synonymous_variants

def _summarize_subset(df: pd.DataFrame, score_col: str) -> Optional[Dict[str, Any]]:
    if score_col not in df.columns or df.empty:
        return None
    scores = pd.to_numeric(df[score_col], errors='coerce').dropna()
    if scores.empty:
        return None
    std = float(scores.std(ddof=1)) if len(scores) > 1 else 0
    return {
        'avg': round(float(scores.mean()), 3),
        'median': round(float(scores.median()), 3),
        'min': round(float(scores.min()), 3),
        'max': round(float(scores.max()), 3),
        'std': round(std, 2),
    }


# TODO: #47 Separate DNA-level and AA-level stats into different functions or modes
def calculate_summary_stats(scores_df: pd.DataFrame, score_col: str) -> Dict[str, Any]:
    """
    Calculate summary stats from variant scores.
    
    This function computes stats for different variant categories including
    overall stats, wild-type scores, synonymous variants, nonsense variants,
    and missense variants.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing variant scores and annotations
    score_col : str
        Score column name to use.
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing calculated summary stats
        
    Examples
    --------
    >>> stats = calculate_summary_stats(scores_df, 'score')
    >>> print(stats['overall'])
    """
    stats: Dict[str, Any] = {}

    if score_col not in scores_df.columns:
        logging.warning(f"Score column '{score_col}' not found in data")
        return stats

    overall_summary = _summarize_subset(scores_df, score_col)
    if overall_summary is not None:
        stats['overall'] = overall_summary

    # Add WT-like reference stats if available.
    if 'annotate_dna' in scores_df.columns:
        wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
        wt_summary = _summarize_subset(wt_subset, score_col)

        synonymous_wt_subset = scores_df[scores_df['annotate_dna'] == 'synonymous']
        synonymous_wt_summary = _summarize_subset(synonymous_wt_subset, score_col)

        if wt_summary is not None:
            stats['wt'] = wt_summary
        if synonymous_wt_summary is not None:
            stats['synonymous_wt'] = synonymous_wt_summary
    elif 'annotate_aa' in scores_df.columns:
        synonymous_wt_subset = scores_df[scores_df['annotate_aa'] == 'synonymous']
        synonymous_wt_summary = _summarize_subset(synonymous_wt_subset, score_col)
        if synonymous_wt_summary is not None:
            stats['synonymous_wt'] = synonymous_wt_summary

    # Summarize AA-level categories from the file being described.
    aa_scores: Optional[pd.DataFrame] = None
    if 'dna_seq_diff' in scores_df.columns and 'aa_seq_diff' in scores_df.columns:
        aa_scores = aggregate_synonymous_variants(scores_df)
    elif 'annotate_aa' in scores_df.columns:
        aa_scores = scores_df

    if aa_scores is not None and score_col in aa_scores.columns:
        nonsense_subset = aa_scores[aa_scores['annotate_aa'] == 'nonsense']
        nonsense_summary = _summarize_subset(nonsense_subset, score_col)

        missense_subset = aa_scores[aa_scores['annotate_aa'] == 'missense_aa']
        missense_summary = _summarize_subset(missense_subset, score_col)

        if nonsense_summary is not None:
            stats['nonsense'] = nonsense_summary
        if missense_summary is not None:
            stats['missense'] = missense_summary

    return stats


def save_summary_stats(
    stats: Dict[str, Any],
    experiment,
    scores_dir: str,
    analysis_logger,
    *,
    stats_basename: str = "stats",
    output_field: str,
    include_in_log_metadata: bool = False,
) -> str:
    """
    Save summary stats to JSON file with logging.
    
    Parameters
    ----------
    stats : Dict[str, Any]
        Dictionary containing summary stats
    experiment : ExperimentConfig
        Experiment configuration
    scores_dir : str
        Directory to save stats file
    analysis_logger : AnalysisLogger
        Logger instance for recording outputs
    stats_basename : str, optional
        Basename used for the stats JSON filename (default: ``"stats"``).
    output_field : str
        Attribute name on ``analysis_logger.outputs`` to store this output (e.g.
        ``"dna_statistics"`` or ``"aa_statistics"``).
    include_in_log_metadata : bool, optional
        If True, embed the stats dictionary into the analysis log metadata for this
        output. Default False to keep logs small (the JSON file is still written).
        
    Returns
    -------
    str
        Path to saved stats file
        
    Examples
    --------
    >>> stats = calculate_summary_stats(scores_df, 'score')
    >>> stats_file = save_summary_stats(stats, experiment, 'output/scores', logger)
    """
    stats_file = os.path.join(scores_dir, f"{experiment.experiment_name}_{stats_basename}.json")
    
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logging.info(f"Saved stats to {stats_file}")
    
    # Log file output
    metadata = {}
    if include_in_log_metadata:
        metadata["stats_summary"] = stats

    analysis_logger.log_output_file(
        output_field,
        f"{experiment.experiment_name}_{stats_basename}.json",
        stats_file,
        **metadata,
    )
    
    return stats_file
