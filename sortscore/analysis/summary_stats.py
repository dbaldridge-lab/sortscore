"""
Summary stats calculation and export for Sort-seq analysis.

This module provides functions for calculating comprehensive summary stats
from Sort-seq variant analysis results and saving them to JSON format.
"""
import json
import os
import logging
import pandas as pd
from typing import Dict, Any, Optional
from sortscore.analysis.variant_aggregation import aggregate_synonymous_variants

# TODO: #47 Separate DNA-level and AA-level stats into different functions or modes
def calculate_summary_stats(scores_df: pd.DataFrame, experiment, score_col: Optional[str] = None) -> Dict[str, Any]:
    """
    Calculate comprehensive summary stats from variant scores.
    
    This function computes stats for different variant categories including
    overall stats, wild-type scores, synonymous variants, nonsense variants,
    and missense variants.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing variant scores and annotations
    experiment : ExperimentConfig
        Experiment configuration containing metadata
    score_col : Optional[str]
        Score column name to use. If None, determined from experiment.avg_method
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing calculated summary stats
        
    Examples
    --------
    >>> stats = calculate_summary_stats(scores_df, experiment)
    >>> print(stats['all_avg'])  # Overall average score
    """
    stats = {}
    
    # Determine score column if not provided
    if score_col is None:
        if experiment.avg_method == 'simple-avg':
            score_col = 'avgscore'
        else:
            score_col_suffix = experiment.avg_method.replace('-', '_')
            score_col = f'avgscore_{score_col_suffix}'
    
    if score_col not in scores_df.columns:
        logging.warning(f"Score column '{score_col}' not found in data")
        return stats
    
    # Overall stats
    mean_val = scores_df[score_col].mean()
    min_val = scores_df[score_col].min()
    max_val = scores_df[score_col].max()
    stats['all_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
    stats['all_min'] = round(float(min_val)) if pd.notna(min_val) else None
    stats['all_max'] = round(float(max_val)) if pd.notna(max_val) else None
    
    # Add annotation-based stats if available
    if 'annotate_dna' in scores_df.columns:
        # WT stats from DNA level
        wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
        if len(wt_subset) > 0:
            stats['wt_dna'] = round(float(mean_val)) if pd.notna(mean_val) else None
        
        # Synonymous (WT) stats from DNA level
        syn_subset = scores_df[scores_df['annotate_dna'] == 'synonymous']
        if len(syn_subset) > 0:
            mean_val = syn_subset[score_col].mean()
            min_val = syn_subset[score_col].min()
            max_val = syn_subset[score_col].max()
            stats['synonymous_wt_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
            stats['synonymous_wt_min'] = round(float(min_val)) if pd.notna(min_val) else None
            stats['synonymous_wt_max'] = round(float(max_val)) if pd.notna(max_val) else None
    
    # Missense and nonsense: use AA-aggregated scores with annotate_aa
    if 'aa_seq_diff' in scores_df.columns:
        aa_scores = aggregate_synonymous_variants(scores_df)
        if score_col in aa_scores.columns:
            # Synonymous stats from AA level (aggregated synonymous variants)
            syn_aa_subset = aa_scores[aa_scores['annotate_aa'] == 'synonymous']
            if len(syn_aa_subset) > 0:
                mean_val = syn_aa_subset[score_col].mean()
                min_val = syn_aa_subset[score_col].min()
                max_val = syn_aa_subset[score_col].max()
                #TODO: #47 This currently gets added to DNA stats
                stats['syn_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                stats['syn_min'] = round(float(min_val)) if pd.notna(min_val) else None
                stats['syn_max'] = round(float(max_val)) if pd.notna(max_val) else None
            
            # Nonsense stats from AA level
            nonsense_subset = aa_scores[aa_scores['annotate_aa'] == 'nonsense']
            if len(nonsense_subset) > 0:
                mean_val = nonsense_subset[score_col].mean()
                min_val = nonsense_subset[score_col].min()
                max_val = nonsense_subset[score_col].max()
                stats['nonsense_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                stats['nonsense_min'] = round(float(min_val)) if pd.notna(min_val) else None
                stats['nonsense_max'] = round(float(max_val)) if pd.notna(max_val) else None
            
            # Missense stats from AA level
            missense_subset = aa_scores[aa_scores['annotate_aa'] == 'missense_aa']
            if len(missense_subset) > 0:
                mean_val = missense_subset[score_col].mean()
                min_val = missense_subset[score_col].min()
                max_val = missense_subset[score_col].max()
                stats['missense_avg'] = round(float(mean_val)) if pd.notna(mean_val) else None
                stats['missense_min'] = round(float(min_val)) if pd.notna(min_val) else None
                stats['missense_max'] = round(float(max_val)) if pd.notna(max_val) else None
    
    return stats


def save_summary_stats(
    stats: Dict[str, Any],
    experiment,
    scores_dir: str,
    output_suffix: str,
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
    output_suffix : str
        Suffix for output filename
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
    >>> stats = calculate_summary_stats(scores_df, experiment)
    >>> stats_file = save_summary_stats(stats, experiment, 'output/scores', 'suffix', logger)
    """
    stats_file = os.path.join(scores_dir, f"{experiment.experiment_name}_{stats_basename}_{output_suffix}.json")
    
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logging.info(f"Saved stats to {stats_file}")
    
    # Log file output
    metadata = {}
    if include_in_log_metadata:
        metadata["stats_summary"] = stats

    analysis_logger.log_output_file(
        output_field,
        f"{experiment.experiment_name}_{stats_basename}_{output_suffix}.json",
        stats_file,
        **metadata,
    )
    
    return stats_file
