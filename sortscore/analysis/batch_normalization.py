"""
Batch normalization module for combining Sort-seq experiments.

This module combines multiple experiments to enable cross-experiment comparisons. 

It supports two normalization approaches:
1. **Z-score scaled 2-pole normalization** (default):
   - Step 1: WT normalization to global reference
   - Step 2: Z-score transformation using synonymous distribution  
   - Step 3: Pathogenic control normalization
   Creates standardized scale where synonymous variants center around 0 with unit 
   variance, making cross-experiment comparisons more easily interpretable.

2. **2-pole normalization**: Uses synonymous and pathogenic variants as reference 
   points with formula: (b/(a-c))*(A-C)
"""

import logging
import os
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from sortscore.utils.load_experiment import ExperimentConfig
from sortscore.analysis.score import calculate_full_activity_scores
from sortscore.analysis.annotation import annotate_scores_dataframe

logger = logging.getLogger(__name__)


def run_batch_analysis(batch_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run complete batch normalization workflow.
    
    Parameters
    ----------
    batch_config : Dict[str, Any]
        Batch configuration containing experiment configs, normalization method, 
        and output settings
        
    Returns
    -------
    Dict[str, Any]
        Results dictionary with combined scores, stats, and metadata
    """
    logger.info("Starting batch normalization workflow")
    
    # Extract configuration parameters
    experiment_configs = batch_config['experiment_configs']
    method = batch_config.get('batch_normalization_method', 'zscore_2pole')
    pathogenic_control_type = batch_config.get('pathogenic_control_type', 'nonsense')
    pathogenic_variants = batch_config.get('pathogenic_variants', None)
    output_dir = batch_config.get('combined_output_dir', '.')
    
    # Load all experiment configurations
    experiments = []
    for config_path in experiment_configs:
        try:
            experiment = ExperimentConfig.from_json(config_path)
            experiments.append(experiment)
            logger.info(f"Loaded experiment config: {config_path}")
        except Exception as e:
            logger.error(f"Failed to load experiment config {config_path}: {e}")
            raise
    
    # Run individual analyses if needed
    all_scores = {}
    all_stats = {}
    
    for i, experiment in enumerate(experiments, 1):
        logger.info(f"Processing experiment {i}: {experiment.experiment_name}")
        
        # Load counts. Avoid re-reading files.
        if experiment.counts is None:
            experiment.load_counts()
        
        # Calculate activity scores
        merged_df = experiment.get_merged_counts()
        scores_df = calculate_full_activity_scores(
            counts=experiment.counts,
            mfi=experiment.mfi,
            min_bins=experiment.bins_required,
            min_reps=experiment.reps_required,
            minread_threshold=experiment.minread_threshold,
            avg_method=experiment.avg_method,
            total_reads=experiment.total_reads,
            cell_prop=experiment.cell_prop,
            merged_df=merged_df,
            max_cv=experiment.max_cv
        )
        
        # Annotate sequences
        scores_df = annotate_scores_dataframe(
            scores_df, experiment.wt_seq, experiment.variant_type
        )
        
        all_scores[f'experiment{i}'] = scores_df
        
        # Extract key statistics
        stats = extract_experiment_stats(scores_df, experiment.avg_method)
        all_stats[f'experiment{i}'] = stats
    
    # Combine raw scores with batch tracking
    combined_scores = combine_raw_scores(all_scores)
    
    # Apply selected normalization method
    if method == '2pole':
        normalized_scores, combined_stats = apply_2pole_normalization(
            combined_scores, all_stats, experiments, pathogenic_control_type, pathogenic_variants
        )
    elif method == 'zscore_2pole':
        normalized_scores, combined_stats = apply_zscore_2pole_normalization(
            combined_scores, all_stats, experiments, pathogenic_control_type, pathogenic_variants
        )
    elif method == 'zscore_center':
        normalized_scores, combined_stats = apply_zscore_center_normalization(
            combined_scores, all_stats, experiments
        )
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    logger.info(f"Applied {method} normalization to {len(normalized_scores)} variants")
    
    return {
        'normalized_scores': normalized_scores,
        'combined_stats': combined_stats,
        'method': method,
        'experiments': experiments,
        'output_dir': output_dir
    }


def combine_raw_scores(all_scores: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Combine raw scores from all experiments with batch tracking.
    
    Parameters
    ----------
    all_scores : Dict[str, pd.DataFrame]
        Dictionary mapping experiment names to score DataFrames
        
    Returns
    -------
    pd.DataFrame
        Combined scores DataFrame with batch column
    """
    combined_dfs = []
    
    for batch_name, scores_df in all_scores.items():
        # Add batch tracking column
        scores_df = scores_df.copy()
        scores_df['batch'] = batch_name
        combined_dfs.append(scores_df)
    
    combined_scores = pd.concat(combined_dfs, ignore_index=True)
    logger.info(f"Combined {len(combined_scores)} variants from {len(all_scores)} experiments")
    
    return combined_scores


def extract_experiment_stats(scores_df: pd.DataFrame, avg_method: str) -> Dict[str, float]:
    """
    Extract key statistics from experiment scores.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        Scores DataFrame for single experiment
    avg_method : str
        Averaging method used for scores
        
    Returns
    -------
    Dict[str, float]
        Dictionary of key statistics
    """
    stats = {}
    score_col_suffix = avg_method.replace('-', '_')
    score_col = f'avgscore_{score_col_suffix}'
    
    if score_col not in scores_df.columns:
        logger.warning(f"Score column {score_col} not found, using 'avgscore'")
        score_col = 'avgscore'
    
    # WT DNA score (if available)
    if 'annotate_dna' in scores_df.columns:
        wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
        if len(wt_subset) > 0:
            stats['wt_dna_score'] = float(wt_subset[score_col].mean())
    
    # Synonymous median
    if 'annotate_dna' in scores_df.columns:
        syn_subset = scores_df[scores_df['annotate_dna'] == 'synonymous']
        if len(syn_subset) > 0:
            stats['syn_median'] = float(syn_subset[score_col].median())
    
    # Nonsense average  
    if 'annotate_aa' in scores_df.columns:
        nonsense_subset = scores_df[scores_df['annotate_aa'] == 'nonsense']
        if len(nonsense_subset) > 0:
            stats['non_avg'] = float(nonsense_subset[score_col].mean())
    
    # Missense average
    if 'annotate_aa' in scores_df.columns:
        missense_subset = scores_df[scores_df['annotate_aa'] == 'missense_aa']
        if len(missense_subset) > 0:
            stats['missense_avg'] = float(missense_subset[score_col].mean())
    
    return stats


def apply_2pole_normalization(
    combined_scores: pd.DataFrame, 
    all_stats: Dict[str, Dict[str, float]],
    experiments: List[Any],
    pathogenic_control_type: str,
    pathogenic_variants: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Apply 2-pole normalization method.
    
    Parameters
    ----------
    combined_scores : pd.DataFrame
        Combined raw scores from all experiments
    all_stats : Dict[str, Dict[str, float]]
        Statistics from individual experiments
    experiments : List[Any]
        List of experiment config objects
    pathogenic_control_type : str
        Type of pathogenic control ('nonsense' or 'custom')
    pathogenic_variants : Optional[List[str]]
        Custom pathogenic variants if not using nonsense
        
    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        Normalized scores DataFrame and combined statistics
    """
    logger.info("Applying 2-pole normalization")
    
    # Calculate global medians
    syn_scores = []
    path_scores = []
    
    for batch_name, batch_df in combined_scores.groupby('batch'):
        # Synonymous scores
        if 'annotate_dna' in batch_df.columns:
            syn_subset = batch_df[batch_df['annotate_dna'] == 'synonymous']
            if len(syn_subset) > 0:
                syn_scores.extend(syn_subset['avgscore'].dropna().tolist())
        
        # Pathogenic scores
        if pathogenic_control_type == 'nonsense':
            if 'annotate_aa' in batch_df.columns:
                path_subset = batch_df[batch_df['annotate_aa'] == 'nonsense']
                if len(path_subset) > 0:
                    path_scores.extend(path_subset['avgscore'].dropna().tolist())
    
    A = np.median(syn_scores) if syn_scores else np.nan
    C = np.median(path_scores) if path_scores else np.nan
    
    # Calculate normalization factors for each experiment
    normalization_factors = {}
    for exp_name, stats in all_stats.items():
        a_i = stats.get('syn_median', np.nan)
        c_i = stats.get('non_avg', np.nan)
        
        if not (np.isnan(a_i) or np.isnan(c_i) or np.isnan(A) or np.isnan(C)):
            if (a_i - c_i) != 0:
                normalization_factors[exp_name] = (A - C) / (a_i - c_i)
            else:
                logger.warning(f"Zero denominator for {exp_name}, skipping normalization")
                normalization_factors[exp_name] = 1.0
        else:
            logger.warning(f"Missing values for {exp_name}, using factor 1.0")
            normalization_factors[exp_name] = 1.0
    
    # Apply normalization to scores
    normalized_scores = combined_scores.copy()
    score_columns = [col for col in normalized_scores.columns if 'score' in col.lower()]
    
    for batch_name, batch_df in normalized_scores.groupby('batch'):
        if batch_name in normalization_factors:
            factor = normalization_factors[batch_name]
            batch_mask = normalized_scores['batch'] == batch_name
            
            for col in score_columns:
                if col in normalized_scores.columns:
                    normalized_scores.loc[batch_mask, col] = (
                        normalized_scores.loc[batch_mask, col] * factor
                    )
    
    # Create combined statistics
    combined_stats = {
        'syn_median_global_raw': A,
        'non_median_global_raw': C,
        'normalization_method': '2pole',
        'pathogenic_control_type': pathogenic_control_type
    }
    
    # Add experiment-specific stats and factors
    for exp_name, stats in all_stats.items():
        combined_stats[f'syn_median_{exp_name}_raw'] = stats.get('syn_median', np.nan)
        combined_stats[f'non_avg_{exp_name}_raw'] = stats.get('non_avg', np.nan)
        combined_stats[f'normalization_factor_{exp_name}'] = normalization_factors.get(exp_name, 1.0)
    
    # Recalculate final statistics from normalized data
    final_stats = recalculate_final_statistics(normalized_scores, 'avgscore')
    combined_stats.update(final_stats)
    
    logger.info(f"2-pole normalization complete with factors: {normalization_factors}")
    return normalized_scores, combined_stats


def apply_zscore_2pole_normalization(
    combined_scores: pd.DataFrame,
    all_stats: Dict[str, Dict[str, float]],
    experiments: List[Any],
    pathogenic_control_type: str,
    pathogenic_variants: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Apply z-score scaled 2-pole normalization method.
    
    Parameters
    ----------
    combined_scores : pd.DataFrame
        Combined raw scores from all experiments
    all_stats : Dict[str, Dict[str, float]]
        Statistics from individual experiments
    experiments : List[Any]
        List of experiment config objects
    pathogenic_control_type : str
        Type of pathogenic control ('nonsense' or 'custom')
    pathogenic_variants : Optional[List[str]]
        Custom pathogenic variants if not using nonsense
        
    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        Normalized scores DataFrame and combined statistics
    """
    logger.info("Applying z-score scaled 2-pole normalization")
    
    # Calculate global WT average
    wt_scores = []
    for batch_name, batch_df in combined_scores.groupby('batch'):
        if 'annotate_dna' in batch_df.columns:
            wt_subset = batch_df[batch_df['annotate_dna'] == 'wt_dna']
            if len(wt_subset) > 0:
                wt_scores.extend(wt_subset['avgscore'].dropna().tolist())
    
    A = np.mean(wt_scores) if wt_scores else np.nan
    
    # Calculate WT normalization factors
    wt_normalization_factors = {}
    for exp_name, stats in all_stats.items():
        a_i = stats.get('wt_dna_score', np.nan)
        if not np.isnan(a_i) and a_i != 0 and not np.isnan(A):
            wt_normalization_factors[exp_name] = A / a_i
        else:
            logger.warning(f"Invalid WT score for {exp_name}, using factor 1.0")
            wt_normalization_factors[exp_name] = 1.0
    
    # Apply WT normalization
    wt_normalized_scores = combined_scores.copy()
    score_columns = [col for col in wt_normalized_scores.columns if 'score' in col.lower()]
    
    for batch_name, batch_df in wt_normalized_scores.groupby('batch'):
        if batch_name in wt_normalization_factors:
            factor = wt_normalization_factors[batch_name]
            batch_mask = wt_normalized_scores['batch'] == batch_name
            
            for col in score_columns:
                if col in wt_normalized_scores.columns:
                    wt_normalized_scores.loc[batch_mask, col] = (
                        wt_normalized_scores.loc[batch_mask, col] * factor
                    )
    
    # Calculate synonymous distribution and apply z-score transformation
    syn_scores = []
    for batch_name, batch_df in wt_normalized_scores.groupby('batch'):
        if 'annotate_dna' in batch_df.columns:
            syn_subset = batch_df[batch_df['annotate_dna'] == 'synonymous']
            if len(syn_subset) > 0:
                syn_scores.extend(syn_subset['avgscore'].dropna().tolist())
    
    syn_mean = np.mean(syn_scores) if syn_scores else 0.0
    syn_std = np.std(syn_scores, ddof=1) if len(syn_scores) > 1 else 1.0
    
    # Apply Z-score transformation
    zscore_normalized_scores = wt_normalized_scores.copy()
    for col in score_columns:
        if col in zscore_normalized_scores.columns and syn_std != 0:
            zscore_normalized_scores[col] = (zscore_normalized_scores[col] - syn_mean) / syn_std
    
    # Recalculate pathogenic controls from z-score normalized values
    zscore_stats = {}
    for batch_name, batch_df in zscore_normalized_scores.groupby('batch'):
        stats = {}
        
        # Pathogenic control average
        if pathogenic_control_type == 'nonsense':
            if 'annotate_aa' in batch_df.columns:
                path_subset = batch_df[batch_df['annotate_aa'] == 'nonsense']
                if len(path_subset) > 0:
                    stats['non_avg_zscore'] = float(path_subset['avgscore'].mean())
        
        zscore_stats[batch_name] = stats
    
    # Calculate global pathogenic average from z-score normalized values
    path_zscore_scores = []
    for stats in zscore_stats.values():
        if 'non_avg_zscore' in stats:
            path_zscore_scores.append(stats['non_avg_zscore'])
    
    C = np.mean(path_zscore_scores) if path_zscore_scores else 0.0
    
    # Calculate pathogenic normalization factors
    path_normalization_factors = {}
    for exp_name, stats in zscore_stats.items():
        c_i = stats.get('non_avg_zscore', np.nan)
        if not np.isnan(c_i) and c_i != 0:
            path_normalization_factors[exp_name] = C / c_i
        else:
            logger.warning(f"Invalid pathogenic score for {exp_name}, using factor 1.0")
            path_normalization_factors[exp_name] = 1.0
    
    # Apply pathogenic normalization
    final_scores = zscore_normalized_scores.copy()
    for batch_name, batch_df in final_scores.groupby('batch'):
        if batch_name in path_normalization_factors:
            factor = path_normalization_factors[batch_name]
            batch_mask = final_scores['batch'] == batch_name
            
            for col in score_columns:
                if col in final_scores.columns:
                    final_scores.loc[batch_mask, col] = (
                        final_scores.loc[batch_mask, col] * factor
                    )
    
    # Create combined statistics
    combined_stats = {
        'wt_dna_score_global_raw': A,
        'syn_mean_global_zscore': syn_mean,
        'syn_std_global_zscore': syn_std,
        'non_avg_global_zscore_z': C,
        'normalization_method': 'zscore_2pole',
        'pathogenic_control_type': pathogenic_control_type
    }
    
    # Add experiment-specific stats and factors
    for exp_name, stats in all_stats.items():
        combined_stats[f'wt_dna_score_{exp_name}_raw'] = stats.get('wt_dna_score', np.nan)
        combined_stats[f'normalization_factor_wt_{exp_name}'] = wt_normalization_factors.get(exp_name, 1.0)
    
    for exp_name, stats in zscore_stats.items():
        combined_stats[f'non_avg_{exp_name}_zscore_z'] = stats.get('non_avg_zscore', np.nan)
        combined_stats[f'normalization_factor_path_{exp_name}'] = path_normalization_factors.get(exp_name, 1.0)
    
    # Recalculate final statistics from normalized data
    final_stats = recalculate_final_statistics(final_scores, 'avgscore')
    combined_stats.update(final_stats)
    
    logger.info("Z-score 2-pole normalization complete")
    return final_scores, combined_stats


def apply_zscore_center_normalization(
    combined_scores: pd.DataFrame,
    all_stats: Dict[str, Dict[str, float]],
    experiments: List[Any]
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Apply z-score centering normalization method using only WT/synonymous variants.
    
    This method is useful when pathogenic controls are not available. It performs
    WT normalization followed by z-score transformation using synonymous variants.
    For DNA variants, uses wt_dna scores; for AA variants, uses synonymous variants.
    
    Parameters
    ----------
    combined_scores : pd.DataFrame
        Combined raw scores from all experiments
    all_stats : Dict[str, Dict[str, float]]
        Statistics from individual experiments
    experiments : List[Any]
        List of experiment config objects to determine variant types
        
    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        Normalized scores DataFrame and combined statistics
    """
    logger.info("Applying z-score centering normalization (WT-only)")
    
    # Determine variant type from experiment configs
    variant_types = [exp.variant_type for exp in experiments]
    if len(set(variant_types)) > 1:
        raise ValueError(f"All experiments must have the same variant_type. Found: {set(variant_types)}")
    variant_type = variant_types[0]
    
    # Calculate global reference scores for normalization
    reference_scores = []
    reference_type = None
    
    if variant_type in {'codon', 'snv', 'dna'}:
        # For DNA variants, prefer wt_dna if available, otherwise use synonymous
        for batch_name, batch_df in combined_scores.groupby('batch'):
            if 'annotate_dna' in batch_df.columns:
                wt_subset = batch_df[batch_df['annotate_dna'] == 'wt_dna']
                if len(wt_subset) > 0:
                    reference_scores.extend(wt_subset['avgscore'].dropna().tolist())
                    reference_type = 'wt_dna'
        
        # Fallback to synonymous if no wt_dna found
        if not reference_scores:
            for batch_name, batch_df in combined_scores.groupby('batch'):
                if 'annotate_dna' in batch_df.columns:
                    syn_subset = batch_df[batch_df['annotate_dna'] == 'synonymous']
                    if len(syn_subset) > 0:
                        reference_scores.extend(syn_subset['avgscore'].dropna().tolist())
                        reference_type = 'synonymous'
    else:
        # For AA variants, use synonymous variants as reference
        for batch_name, batch_df in combined_scores.groupby('batch'):
            if 'annotate_dna' in batch_df.columns:
                syn_subset = batch_df[batch_df['annotate_dna'] == 'synonymous']
                if len(syn_subset) > 0:
                    reference_scores.extend(syn_subset['avgscore'].dropna().tolist())
                    reference_type = 'synonymous'
    
    if not reference_scores:
        raise ValueError("No reference variants found for normalization. Need either wt_dna or synonymous variants.")
    
    global_reference_avg = np.mean(reference_scores)
    logger.info(f"Using {reference_type} variants as reference with global mean: {global_reference_avg:.3f}")
    
    # Calculate normalization factors based on variant type
    normalization_factors = {}
    for exp_name, stats in all_stats.items():
        if reference_type == 'wt_dna':
            exp_reference_score = stats.get('wt_dna_score', np.nan)
        else:  # synonymous
            exp_reference_score = stats.get('syn_median', np.nan)
        
        if not np.isnan(exp_reference_score) and exp_reference_score != 0 and not np.isnan(global_reference_avg):
            normalization_factors[exp_name] = global_reference_avg / exp_reference_score
        else:
            logger.warning(f"Invalid {reference_type} score for {exp_name}, using factor 1.0")
            normalization_factors[exp_name] = 1.0
    
    # Apply reference normalization
    reference_normalized_scores = combined_scores.copy()
    score_columns = [col for col in reference_normalized_scores.columns if 'score' in col.lower()]
    
    for batch_name, batch_df in reference_normalized_scores.groupby('batch'):
        if batch_name in normalization_factors:
            factor = normalization_factors[batch_name]
            batch_mask = reference_normalized_scores['batch'] == batch_name
            
            for col in score_columns:
                if col in reference_normalized_scores.columns:
                    reference_normalized_scores.loc[batch_mask, col] = (
                        reference_normalized_scores.loc[batch_mask, col] * factor
                    )
    
    # Calculate synonymous distribution for z-score transformation
    syn_scores = []
    for batch_name, batch_df in reference_normalized_scores.groupby('batch'):
        if 'annotate_dna' in batch_df.columns:
            syn_subset = batch_df[batch_df['annotate_dna'] == 'synonymous']
            if len(syn_subset) > 0:
                syn_scores.extend(syn_subset['avgscore'].dropna().tolist())
    
    if not syn_scores:
        logger.warning("No synonymous variants found for z-score transformation, using all variants")
        syn_scores = reference_normalized_scores['avgscore'].dropna().tolist()
    
    syn_mean = np.mean(syn_scores) if syn_scores else 0.0
    syn_std = np.std(syn_scores, ddof=1) if len(syn_scores) > 1 else 1.0
    
    # Apply z-score transformation (final normalization)
    final_scores = reference_normalized_scores.copy()
    for col in score_columns:
        if col in final_scores.columns and syn_std != 0:
            final_scores[col] = (final_scores[col] - syn_mean) / syn_std
    
    # Create combined statistics (no pathogenic control type since not used)
    combined_stats = {
        f'{reference_type}_score_global_raw': global_reference_avg,
        'syn_mean_global_ref_norm': syn_mean,
        'syn_std_global_ref_norm': syn_std,
        'normalization_method': 'zscore_center',
        'reference_type': reference_type
    }
    
    # Add experiment-specific stats and factors
    for exp_name, stats in all_stats.items():
        if reference_type == 'wt_dna':
            combined_stats[f'wt_dna_score_{exp_name}_raw'] = stats.get('wt_dna_score', np.nan)
        else:
            combined_stats[f'syn_median_{exp_name}_raw'] = stats.get('syn_median', np.nan)
        combined_stats[f'normalization_factor_{exp_name}'] = normalization_factors.get(exp_name, 1.0)
    
    # Recalculate final statistics from normalized data
    final_stats = recalculate_final_statistics(final_scores, 'avgscore')
    combined_stats.update(final_stats)
    
    logger.info("Z-score centering normalization complete")
    return final_scores, combined_stats


def recalculate_final_statistics(
    normalized_scores: pd.DataFrame, 
    score_col: str
) -> Dict[str, Any]:
    """
    Recalculate statistics from fully normalized score data.
    
    Parameters
    ----------
    normalized_scores : pd.DataFrame
        Fully normalized scores DataFrame
    score_col : str
        Name of main score column to use for calculations
        
    Returns
    -------
    Dict[str, Any]
        Dictionary of recalculated statistics
    """
    final_stats = {}
    
    # Global final statistics
    if score_col in normalized_scores.columns:
        scores = normalized_scores[score_col].dropna()
        if len(scores) > 0:
            final_stats['all_avg_global_final'] = float(scores.mean())
            final_stats['all_min_global_final'] = float(scores.min())
            final_stats['all_max_global_final'] = float(scores.max())
    
    # Per-batch final statistics
    for batch_name, batch_df in normalized_scores.groupby('batch'):
        batch_prefix = f'{batch_name}_final'
        
        if score_col in batch_df.columns:
            scores = batch_df[score_col].dropna()
            if len(scores) > 0:
                final_stats[f'all_avg_{batch_prefix}'] = float(scores.mean())
        
        # Missense average per batch
        if 'annotate_aa' in batch_df.columns:
            missense_subset = batch_df[batch_df['annotate_aa'] == 'missense_aa']
            if len(missense_subset) > 0 and score_col in missense_subset.columns:
                missense_scores = missense_subset[score_col].dropna()
                if len(missense_scores) > 0:
                    final_stats[f'missense_avg_{batch_prefix}'] = float(missense_scores.mean())
        
        # Nonsense average per batch  
        if 'annotate_aa' in batch_df.columns:
            nonsense_subset = batch_df[batch_df['annotate_aa'] == 'nonsense']
            if len(nonsense_subset) > 0 and score_col in nonsense_subset.columns:
                nonsense_scores = nonsense_subset[score_col].dropna()
                if len(nonsense_scores) > 0:
                    final_stats[f'nonsense_avg_{batch_prefix}'] = float(nonsense_scores.mean())
    
    return final_stats


def generate_batch_visualizations(
    results: Dict[str, Any],
    batch_config: Dict[str, Any],
    output_dir: str,
    suffix: Optional[str] = None
) -> None:
    """
    Generate tiled heatmap visualizations from batch results.
    
    Parameters
    ----------
    results : Dict[str, Any]
        Results dictionary from batch analysis
    batch_config : Dict[str, Any] 
        Batch configuration dictionary
    output_dir : str
        Output directory for visualizations
    suffix : Optional[str]
        Optional suffix for output files
    """
    from sortscore.analysis.batch_config import BatchConfig
    from sortscore.visualization.heatmaps import plot_tiled_heatmap
    from sortscore.utils.file_utils import ensure_output_subdirs
    
    # Ensure figures directory exists
    ensure_output_subdirs(output_dir)
    figures_dir = os.path.join(output_dir, 'figures')
    
    # Create batch config object for visualization
    config_obj = BatchConfig(**batch_config)
    
    # Generate suffix if not provided
    if suffix is None:
        timestamp = pd.Timestamp.now().strftime('%Y%m%d')
        method = results['method']
        suffix = f"batch_{method}_{timestamp}"
    
    # Determine score column to visualize
    score_col = 'avgscore'  # Default, could be made configurable
    
    # Calculate WT score for colorbar reference if available
    wt_score = None
    combined_stats = results.get('combined_stats', {})
    
    # Try to get normalized WT score from stats
    if results['method'] == 'zscore_2pole':
        # For z-score method, WT should be close to 0 after normalization
        wt_score = 0.0
    elif results['method'] == '2pole':
        # For 2-pole method, try to get global synonymous median
        wt_score = combined_stats.get('syn_median_global_raw')
    
    # Create tiled heatmap
    heatmap_file = os.path.join(figures_dir, f"tiled_heatmap_{suffix}.png")
    
    try:
        plot_tiled_heatmap(
            batch_data=results['normalized_scores'],
            score_col=score_col,
            batch_config=config_obj,
            experiments=results['experiments'],
            wt_score=wt_score,
            export=True,
            output=heatmap_file,
            export_matrix=True,
            show_biophysical_properties=True  # Enable biophysical properties panel
        )
        logger.info(f"Generated tiled heatmap: {heatmap_file}")
        
    except Exception as e:
        logger.error(f"Failed to generate tiled heatmap: {e}")


def save_batch_results(
    results: Dict[str, Any],
    output_dir: str,
    suffix: Optional[str] = None
) -> None:
    """
    Save batch normalization results to files.
    
    Parameters
    ----------
    results : Dict[str, Any]
        Results dictionary from batch analysis
    output_dir : str
        Output directory for saving results
    suffix : Optional[str]
        Optional suffix for output files
    """
    from sortscore.utils.file_utils import ensure_output_subdirs
    
    # Ensure output directories exist
    ensure_output_subdirs(output_dir)
    
    # Generate suffix if not provided
    if suffix is None:
        timestamp = pd.Timestamp.now().strftime('%Y%m%d')
        method = results['method']
        suffix = f"batch_{method}_{timestamp}"
    
    # Save normalized scores
    scores_dir = os.path.join(output_dir, 'scores')
    scores_file = os.path.join(scores_dir, f"batch_scores_{suffix}.csv")
    
    normalized_scores = results['normalized_scores'].copy()
    # Round score columns to integers for consistency
    score_columns = [col for col in normalized_scores.columns if 'score' in col.lower()]
    for col in score_columns:
        if normalized_scores[col].dtype in ['float64', 'float32']:
            normalized_scores[col] = normalized_scores[col].round().astype('Int64')
    
    normalized_scores.to_csv(scores_file, index=False)
    logger.info(f"Saved batch scores to {scores_file}")
    
    # Save combined statistics
    stats_file = os.path.join(scores_dir, f"batch_stats_{suffix}.json")
    with open(stats_file, 'w') as f:
        json.dump(results['combined_stats'], f, indent=2)
    logger.info(f"Saved batch statistics to {stats_file}")
