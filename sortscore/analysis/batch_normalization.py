"""
Batch normalization module for combining Sort-seq experiments.

This module combines multiple experiments to enable cross-experiment comparisons. 

It supports two normalization approaches:
1. **Z-score scaled 2-pole normalization** (default):
   - Step 1: Synonymous-variant normalization to global reference
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
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, cast

logger = logging.getLogger(__name__)

from sortscore.analysis.aa_scores import _get_score_column_from_avg_method
from sortscore.analysis.summary_stats import calculate_summary_stats
def _ensure_avgscore_column(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Ensure a canonical `avgscore` column exists for normalization steps."""
    if 'avgscore' in scores_df.columns:
        return scores_df
    for col in ('avgscore_rep_weighted', 'avgscore_simple_avg'):
        if col in scores_df.columns:
            out = scores_df.copy()
            out['avgscore'] = out[col]
            return out
    raise ValueError("Scores file missing required score column ('avgscore' or known avgscore_* variants)")


def _score_columns(scores_df: pd.DataFrame) -> List[str]:
    """Return score-like columns used in normalization math."""
    return [col for col in scores_df.columns if 'score' in col.lower()]


def _coerce_score_columns_to_float(scores_df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure score columns are float-typed before applying multiplicative scaling.

    Pandas warns when assigning float values back into integer-typed columns.
    """
    out = scores_df.copy()
    for col in _score_columns(out):
        out[col] = pd.to_numeric(out[col], errors='coerce').astype(float)
    return out


def _build_stage_section(
    *,
    global_values: Optional[Dict[str, Any]] = None,
    tile_values: Optional[Dict[str, Dict[str, Any]]] = None,
    factor_values: Optional[Dict[str, float]] = None,
    factor_key: str = 'normalization_factor',
) -> Dict[str, Any]:
    stage: Dict[str, Any] = {}
    if global_values:
        stage['global'] = dict(global_values)

    tile_keys = set()
    if tile_values:
        tile_keys.update(tile_values.keys())
    if factor_values:
        tile_keys.update(factor_values.keys())

    for batch_key in sorted(tile_keys):
        tile_stats: Dict[str, Any] = {}
        if tile_values and batch_key in tile_values:
            tile_stats.update(tile_values[batch_key])
        if factor_values and batch_key in factor_values:
            tile_stats[factor_key] = factor_values[batch_key]
        stage[batch_key] = tile_stats

    return stage


def _build_final_batch_stats(normalized_scores: pd.DataFrame) -> Dict[str, Any]:
    """Create final stats for global and per-tile normalized outputs."""
    final_stats: Dict[str, Any] = {
        'global': calculate_summary_stats(normalized_scores, 'avgscore'),
    }
    for batch_name, batch_df in normalized_scores.groupby('batch'):
        final_stats[cast(str, batch_name)] = calculate_summary_stats(
            batch_df,
            'avgscore',
        )
    return final_stats


def _get_heatmap_reference_marker(track_stats: Dict[str, Any]) -> Optional[float]:
    """Return the synonymous reference score for one normalized score track."""
    final_global_stats = track_stats.get('final', {}).get('global', {})
    syn_wt_stats = final_global_stats.get('synonymous_wt')
    if syn_wt_stats is not None and 'avg' in syn_wt_stats:
        return syn_wt_stats['avg']
    return None


def _build_normalization_stats(
    *,
    raw_tile_values: Optional[Dict[str, Dict[str, Any]]] = None,
    raw_global_values: Optional[Dict[str, Any]] = None,
    wt_stage_global_values: Optional[Dict[str, Any]] = None,
    zscore_tile_values: Optional[Dict[str, Dict[str, Any]]] = None,
    zscore_global_values: Optional[Dict[str, Any]] = None,
    final_scores: pd.DataFrame,
    wt_factors: Optional[Dict[str, float]] = None,
    path_factors: Optional[Dict[str, float]] = None,
    normalization_factors: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    stats: Dict[str, Any] = {}

    stats['raw'] = _build_stage_section(
        global_values=raw_global_values,
        tile_values=raw_tile_values,
    )

    if wt_stage_global_values is not None or wt_factors is not None:
        stats['wt_alignment'] = _build_stage_section(
            global_values=wt_stage_global_values,
            factor_values=wt_factors,
        )

    if zscore_global_values is not None or zscore_tile_values is not None or path_factors is not None:
        stats['zscore'] = _build_stage_section(
            global_values=zscore_global_values,
            tile_values=zscore_tile_values,
            factor_values=path_factors,
            factor_key='pathogenic_normalization_factor',
        )

    if normalization_factors is not None:
        scaling_stage = {'per_tile': {}}
        for batch_name, factor in normalization_factors.items():
            scaling_stage['per_tile'][batch_name] = {'normalization_factor': factor}
        stats['scaling'] = scaling_stage

    stats['final'] = _build_final_batch_stats(final_scores)
    return stats


def _load_score_tables_from_output_dir(output_dir: str) -> Dict[str, pd.DataFrame]:
    """Load available DNA and AA score CSVs from `output_dir/scores`."""
    scores_dir = Path(output_dir).expanduser().resolve() / 'scores'
    if not scores_dir.exists():
        raise FileNotFoundError(f"Scores directory not found: {scores_dir}")

    dna_files = sorted(scores_dir.glob("*_dna_scores.csv"))
    aa_files = sorted(scores_dir.glob("*_aa_scores.csv"))
    if not dna_files and not aa_files:
        raise FileNotFoundError(f"No score CSV found in {scores_dir}")

    score_tables: Dict[str, pd.DataFrame] = {}
    if dna_files:
        score_tables['dna'] = _ensure_avgscore_column(pd.read_csv(dna_files[-1]))
    if aa_files:
        score_tables['aa'] = _ensure_avgscore_column(pd.read_csv(aa_files[-1]))
    return score_tables


def _run_normalization(
    all_scores: Dict[str, pd.DataFrame],
    all_stats: Dict[str, Dict[str, float]],
    experiments: List[Any],
    method: str,
    pathogenic_control_type: str,
    pathogenic_variants: Optional[List[str]],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Normalize one set of score tables using the selected method."""
    combined_scores = combine_raw_scores(all_scores)

    if method == 'linear_range':
        return apply_linear_range_normalization(
            combined_scores, all_stats, experiments, pathogenic_control_type, pathogenic_variants
        )
    if method == 'zscore_2pole':
        return apply_zscore_2pole_normalization(
            combined_scores, all_stats, experiments, pathogenic_control_type, pathogenic_variants
        )
    if method == 'zscore_onepole':
        return apply_zscore_onepole_normalization(
            combined_scores, all_stats, experiments
        )
    raise ValueError(f"Unknown normalization method: {method}")


def _get_position_range_from_config(experiment_cfg: Dict[str, Any], experiment_idx: int) -> Tuple[int, int]:
    """Read and validate min/max positions from a tile entry in `batch_config.experiments`."""
    if 'min_pos' not in experiment_cfg or 'max_pos' not in experiment_cfg:
        raise ValueError(
            f"experiments[{experiment_idx}] missing 'min_pos'/'max_pos'; "
            "position range must be set in batch_config tile entries"
        )
    try:
        min_pos = int(experiment_cfg['min_pos'])
        max_pos = int(experiment_cfg['max_pos'])
    except Exception as e:
        raise ValueError(
            f"experiments[{experiment_idx}].min_pos/max_pos must be int-like"
        ) from e
    if min_pos >= max_pos:
        raise ValueError(
            f"experiments[{experiment_idx}] min_pos must be less than max_pos"
        )
    return min_pos, max_pos


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
    
    experiment_entries = batch_config.get('experiments') or []
    method = batch_config.get('batch_normalization_method', 'zscore_2pole')
    pathogenic_control_type = batch_config.get('pathogenic_control_type', 'nonsense')
    pathogenic_variants = batch_config.get('pathogenic_variants', None)
    base_output_dir = Path(str(batch_config.get('combined_output_dir', '.'))).expanduser().resolve()
    
    # Load all tile score tables from batch config entries.
    experiments = []
    all_scores: Dict[str, Dict[str, pd.DataFrame]] = {'dna': {}, 'aa': {}}
    all_stats: Dict[str, Dict[str, Dict[str, float]]] = {'dna': {}, 'aa': {}}

    for idx, cfg in enumerate(experiment_entries, 1):
        try:
            tile = int(cfg['tile'])
            output_dir_i = str(Path(str(cfg['output_dir'])).expanduser().resolve())
            score_tables = _load_score_tables_from_output_dir(output_dir_i)
            batch_name = f"tile{idx}"

            mutagenesis_type = cfg.get('mutagenesis_type')
            if mutagenesis_type is None:
                dna_scores_df = score_tables.get('dna')
                mutagenesis_type = 'codon' if dna_scores_df is not None else 'aa'

            avg_method = cfg.get('avg_method', 'rep-weighted')
            min_pos, max_pos = _get_position_range_from_config(cfg, idx - 1)
            for track_name, scores_df in score_tables.items():
                all_scores[track_name][batch_name] = scores_df
                all_stats[track_name][batch_name] = extract_experiment_stats(scores_df, avg_method)
            experiments.append(
                SimpleNamespace(
                    tile=tile,
                    wt_seq=cfg.get('wt_seq'),
                    mutagenesis_type=mutagenesis_type,
                    avg_method=avg_method,
                    min_pos=min_pos,
                    max_pos=max_pos,
                    mutagenesis_variants=cfg.get('mutagenesis_variants'),
                )
            )
            logger.info(f"Loaded batch config entry {idx}: tile={tile}, output_dir={output_dir_i}")
        except Exception as e:
            logger.error(f"Failed to load batch config entry {idx}: {e}")
            raise

    output_dir = str((base_output_dir / 'normalized' / method).resolve())
    normalized_results: Dict[str, pd.DataFrame] = {'dna': pd.DataFrame(), 'aa': pd.DataFrame()}
    combined_stats: Dict[str, Dict[str, Any]] = {'dna': {}, 'aa': {}}

    for track_name in ('dna', 'aa'):
        if not all_scores[track_name]:
            continue
        normalized_results[track_name], combined_stats[track_name] = _run_normalization(
            all_scores[track_name],
            all_stats[track_name],
            experiments,
            method,
            pathogenic_control_type,
            pathogenic_variants,
        )
        logger.info(
            f"Normalized {len(normalized_results[track_name])} {track_name.upper()} rows with {method}"
        )

    return {
        'normalized_scores': normalized_results['dna'],
        'normalized_aa_scores': normalized_results['aa'],
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
    if 'annotate_aa' not in scores_df.columns:
        raise ValueError("scores_df must contain 'annotate_aa' for synonymous reference selection")
    
    # Synonymous median
    syn_subset = scores_df[scores_df['annotate_aa'] == 'synonymous']
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


def apply_linear_range_normalization(
    combined_scores: pd.DataFrame, 
    all_stats: Dict[str, Dict[str, float]],
    experiments: List[Any],
    pathogenic_control_type: str,
    pathogenic_variants: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Apply linear-range normalization method.
    
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
    logger.info("Applying linear-range normalization")
    if 'annotate_aa' not in combined_scores.columns:
        raise ValueError("scores_df must contain 'annotate_aa' for synonymous reference selection")
    
    # Calculate global medians
    syn_scores = []
    path_scores = []
    
    for batch_name, batch_df in combined_scores.groupby('batch'):
        # Synonymous scores
        syn_subset = batch_df[batch_df['annotate_aa'] == 'synonymous']
        if len(syn_subset) > 0:
            syn_scores.extend(syn_subset['avgscore'].dropna().tolist())
        
        # Pathogenic scores
        if pathogenic_control_type == 'nonsense':
            if 'annotate_aa' in batch_df.columns:
                path_subset = batch_df[batch_df['annotate_aa'] == 'nonsense']
                if len(path_subset) > 0:
                    path_scores.extend(path_subset['avgscore'].dropna().tolist())
    
    if not syn_scores:
        raise ValueError("No synonymous variants found for normalization.")
    if not path_scores:
        raise ValueError("No pathogenic variants found for normalization.")

    A = np.median(syn_scores)
    C = np.median(path_scores)
    
    # Calculate normalization factors for each experiment
    normalization_factors = {}
    for exp_name, stats in all_stats.items():
        a_i = stats.get('syn_median', np.nan)
        c_i = stats.get('non_avg', np.nan)
        
        if np.isnan(a_i) or np.isnan(c_i):
            raise ValueError(f"Missing normalization values for {exp_name}.")
        if (a_i - c_i) == 0:
            raise ValueError(f"Zero normalization denominator for {exp_name}.")
        normalization_factors[exp_name] = (A - C) / (a_i - c_i)
    
    # Apply normalization to scores
    normalized_scores = _coerce_score_columns_to_float(combined_scores)
    score_columns = _score_columns(normalized_scores)
    
    for batch_name, batch_df in normalized_scores.groupby('batch'):
        if batch_name in normalization_factors:
            factor = normalization_factors[batch_name]
            batch_mask = normalized_scores['batch'] == batch_name
            
            for col in score_columns:
                if col in normalized_scores.columns:
                    normalized_scores.loc[batch_mask, col] = (
                        normalized_scores.loc[batch_mask, col] * factor
                    )
    
    combined_stats = _build_normalization_stats(
        raw_tile_values={
            batch_name: {
                'syn_median': stats.get('syn_median'),
                'non_avg': stats.get('non_avg'),
            }
            for batch_name, stats in all_stats.items()
        },
        raw_global_values={
            'syn_median': A,
            'pathogenic_median': C,
        },
        final_scores=normalized_scores,
        normalization_factors=normalization_factors,
    )
    
    logger.info(f"Linear-range normalization complete with factors: {normalization_factors}")
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
    if 'annotate_aa' not in combined_scores.columns:
        raise ValueError("scores_df must contain 'annotate_aa' for synonymous reference selection")
    
    # Use synonymous variants as the alignment reference for every input type.
    reference_scores = []
    for _, batch_df in combined_scores.groupby('batch'):
        syn_subset = batch_df[batch_df['annotate_aa'] == 'synonymous']
        if len(syn_subset) > 0:
            reference_scores.extend(syn_subset['avgscore'].dropna().tolist())

    if not reference_scores:
        raise ValueError("No synonymous variants found for normalization.")

    A = np.median(reference_scores)

    # Calculate initial alignment factors.
    wt_normalization_factors = {}
    for exp_name, stats in all_stats.items():
        a_i = stats.get('syn_median', np.nan)
        if not np.isnan(a_i) and a_i != 0 and not np.isnan(A):
            wt_normalization_factors[exp_name] = A / a_i
        else:
            raise ValueError(f"Invalid synonymous reference score for {exp_name}.")
    
    # Apply initial reference alignment.
    wt_normalized_scores = _coerce_score_columns_to_float(combined_scores)
    score_columns = _score_columns(wt_normalized_scores)
    
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
        syn_subset = batch_df[batch_df['annotate_aa'] == 'synonymous']
        if len(syn_subset) > 0:
            syn_scores.extend(syn_subset['avgscore'].dropna().tolist())
    
    if len(syn_scores) < 2:
        raise ValueError("At least two synonymous scores are required for z-score normalization.")

    syn_mean = np.mean(syn_scores)
    syn_median = float(np.median(syn_scores))
    syn_std = np.std(syn_scores, ddof=1)
    if not np.isfinite(syn_std) or syn_std == 0:
        raise ValueError("Synonymous scores must have a finite, non-zero standard deviation.")
    
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
    
    if not path_zscore_scores:
        raise ValueError("No pathogenic reference scores found for normalization.")
    C = np.mean(path_zscore_scores)
    
    # Calculate pathogenic normalization factors
    path_normalization_factors = {}
    for exp_name, stats in zscore_stats.items():
        c_i = stats.get('non_avg_zscore', np.nan)
        if not np.isnan(c_i) and c_i != 0:
            path_normalization_factors[exp_name] = C / c_i
        else:
            raise ValueError(f"Invalid pathogenic reference score for {exp_name}.")
    
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
    
    combined_stats = _build_normalization_stats(
        raw_tile_values={
            batch_name: {'syn_median': stats.get('syn_median')}
            for batch_name, stats in all_stats.items()
        },
        raw_global_values={
            'reference_type': 'synonymous',
            'reference_median': A,
        },
        wt_stage_global_values={
            'syn_avg': syn_mean,
            'syn_median': syn_median,
            'syn_std': syn_std,
        },
        zscore_tile_values=zscore_stats,
        zscore_global_values={'non_avg_zscore': C},
        final_scores=final_scores,
        wt_factors=wt_normalization_factors,
        path_factors=path_normalization_factors,
    )
    
    logger.info("Z-score 2-pole normalization complete")
    return final_scores, combined_stats


def apply_zscore_onepole_normalization(
    combined_scores: pd.DataFrame,
    all_stats: Dict[str, Dict[str, float]],
    experiments: List[Any]
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Apply z-score centering normalization using synonymous variants.
    
    This method is useful when pathogenic controls are not available. It performs
    Synonymous-reference alignment followed by a z-score transformation using
    synonymous variants. The same reference is used for DNA and AA inputs.
    
    Parameters
    ----------
    combined_scores : pd.DataFrame
        Combined raw scores from all experiments
    all_stats : Dict[str, Dict[str, float]]
        Statistics from individual experiments
    experiments : List[Any]
        List of experiment config objects (retained for API compatibility)
        
    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Any]]
        Normalized scores DataFrame and combined statistics
    """
    logger.info("Applying z-score centering normalization (synonymous reference)")

    if 'annotate_aa' not in combined_scores.columns:
        raise ValueError("scores_df must contain 'annotate_aa' for synonymous reference selection")
    
    # Calculate global reference scores for normalization
    reference_scores = []
    for _, batch_df in combined_scores.groupby('batch'):
        syn_subset = batch_df[batch_df['annotate_aa'] == 'synonymous']
        if len(syn_subset) > 0:
            reference_scores.extend(syn_subset['avgscore'].dropna().tolist())
    
    if not reference_scores:
        raise ValueError("No synonymous variants found for normalization.")
    
    global_reference_median = np.median(reference_scores)
    logger.info(
        "Using synonymous variants as reference with global median: "
        f"{global_reference_median:.3f}"
    )
    
    # Calculate normalization factors based on variant type
    normalization_factors = {}
    for exp_name, stats in all_stats.items():
        exp_reference_score = stats.get('syn_median', np.nan)
        
        if not np.isnan(exp_reference_score) and exp_reference_score != 0 and not np.isnan(global_reference_median):
            normalization_factors[exp_name] = global_reference_median / exp_reference_score
        else:
            raise ValueError(f"Invalid synonymous reference score for {exp_name}.")
    
    # Apply reference normalization
    reference_normalized_scores = _coerce_score_columns_to_float(combined_scores)
    score_columns = _score_columns(reference_normalized_scores)
    
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
        syn_subset = batch_df[batch_df['annotate_aa'] == 'synonymous']
        if len(syn_subset) > 0:
            syn_scores.extend(syn_subset['avgscore'].dropna().tolist())
    
    if len(syn_scores) < 2:
        raise ValueError("At least two synonymous scores are required for z-score normalization.")

    syn_mean = np.mean(syn_scores)
    syn_median = float(np.median(syn_scores))
    syn_std = np.std(syn_scores, ddof=1)
    if not np.isfinite(syn_std) or syn_std == 0:
        raise ValueError("Synonymous scores must have a finite, non-zero standard deviation.")
    
    # Apply z-score transformation (final normalization)
    final_scores = reference_normalized_scores.copy()
    for col in score_columns:
        if col in final_scores.columns and syn_std != 0:
            final_scores[col] = (final_scores[col] - syn_mean) / syn_std
    
    combined_stats = _build_normalization_stats(
        raw_tile_values={
            batch_name: {'syn_median': stats.get('syn_median')}
            for batch_name, stats in all_stats.items()
        },
        raw_global_values={
            'reference_type': 'synonymous',
            'reference_median': global_reference_median,
        },
        wt_stage_global_values={
            'syn_avg': syn_mean,
            'syn_median': syn_median,
            'syn_std': syn_std,
        },
        final_scores=final_scores,
        normalization_factors=normalization_factors,
    )
    
    logger.info("Z-score centering normalization complete")
    return final_scores, combined_stats


def generate_batch_visualizations(
    results: Dict[str, Any],
    batch_config: Dict[str, Any],
    output_dir: str,
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
    """
    from sortscore.analysis.batch_config import BatchConfig
    from sortscore.visualization.heatmaps import plot_tiled_heatmap
    from sortscore.utils.file_utils import ensure_output_subdirs
    
    # Ensure figures directory exists
    ensure_output_subdirs(output_dir)
    figures_dir = os.path.join(output_dir, 'figures')
    
    # Create batch config object for visualization
    config_obj = BatchConfig(**batch_config)
    
    # Determine score column to visualize
    score_col = 'avgscore'  # Default, could be made configurable
    
    combined_stats = results.get('combined_stats', {})
    dna_combined_stats = combined_stats.get('dna', {})
    aa_combined_stats = combined_stats.get('aa', {})

    def _build_colorbar_ticks(
        score_series: pd.Series,
        track_stats: Dict[str, Any],
        reference_score: Optional[float],
    ) -> Tuple[Optional[List[float]], Optional[List[str]]]:
        series = pd.to_numeric(score_series, errors='coerce').dropna()
        if len(series) == 0:
            return None, None

        data_min = float(series.min())
        data_max = float(series.max())

        pathogenic_tick = None
        if config_obj.pathogenic_control_type == 'nonsense':
            pathogenic_tick = (
                track_stats.get('final', {})
                .get('global', {})
                .get('nonsense', {})
                .get('avg')
            )
        if pathogenic_tick is None:
            pathogenic_tick = (
                track_stats.get('zscore', {})
                .get('global', {})
                .get('non_avg_zscore')
            )

        ticks = [(data_min, f'{data_min:.1f}')]
        if reference_score is not None and pd.notna(reference_score):
            ticks.append((float(reference_score), 'Synonymous Avg'))
        if pathogenic_tick is not None and pd.notna(pathogenic_tick):
            ticks.append((float(pathogenic_tick), '(-) Control Avg'))
        ticks.append((data_max, f'{data_max:.1f}'))

        dedup = {}
        for val, label in ticks:
            dedup[round(float(val), 6)] = (float(val), label)
        ordered = sorted(dedup.values(), key=lambda x: x[0])
        return [val for val, _ in ordered], [label for _, label in ordered]

    tick_values = None
    tick_labels = None
    dna_batch_data = results['normalized_scores']
    dna_reference_score = _get_heatmap_reference_marker(dna_combined_stats)
    if score_col in dna_batch_data.columns:
        tick_values, tick_labels = _build_colorbar_ticks(
            dna_batch_data[score_col],
            dna_combined_stats,
            dna_reference_score,
        )

    # Create tiled heatmap
    tiled_heatmap_file = os.path.join(figures_dir, "tiled_heatmap.png")
    combined_aa_heatmap_file = os.path.join(figures_dir, "combined_aa_heatmap.png")

    try:
        if not dna_batch_data.empty:
            plot_tiled_heatmap(
                batch_data=dna_batch_data,
                score_col=score_col,
                batch_config=config_obj,
                experiments=results['experiments'],
                wt_score=dna_reference_score,
                tick_values=tick_values,
                tick_labels=tick_labels,
                export=True,
                output=tiled_heatmap_file,
                export_matrix=True,
                biophysical_properties=False
            )
            logger.info(f"Generated tiled heatmap: {tiled_heatmap_file}")

        # Also export a combined AA heatmap for cross-tile outputs.
        aa_batch_data = results['normalized_aa_scores']
        if not aa_batch_data.empty:
            aa_tick_values = tick_values
            aa_tick_labels = tick_labels
            aa_reference_score = _get_heatmap_reference_marker(aa_combined_stats)
            if score_col in aa_batch_data.columns:
                aa_tick_values, aa_tick_labels = _build_colorbar_ticks(
                    aa_batch_data[score_col],
                    aa_combined_stats,
                    aa_reference_score,
                )
            aa_experiments = [
                SimpleNamespace(
                    tile=exp.tile,
                    wt_seq=exp.wt_seq,
                    mutagenesis_type='aa',
                    avg_method=exp.avg_method,
                    min_pos=exp.min_pos,
                    max_pos=exp.max_pos,
                    mutagenesis_variants=exp.mutagenesis_variants,
                )
                for exp in results['experiments']
            ]

            plot_tiled_heatmap(
                batch_data=aa_batch_data,
                score_col=score_col,
                batch_config=config_obj,
                experiments=aa_experiments,
                wt_score=aa_reference_score,
                tick_values=aa_tick_values,
                tick_labels=aa_tick_labels,
                export=True,
                output=combined_aa_heatmap_file,
                export_matrix=True,
                biophysical_properties=False,
                title='Combined AA Heatmap'
            )
            logger.info(f"Generated combined AA heatmap: {combined_aa_heatmap_file}")

    except Exception as e:
        logger.error(f"Failed to generate tiled/combined AA heatmap: {e}")


def save_batch_results(
    results: Dict[str, Any],
    output_dir: str,
) -> None:
    """
    Save batch normalization results to files.
    
    Parameters
    ----------
    results : Dict[str, Any]
        Results dictionary from batch analysis
    output_dir : str
        Output directory for saving results
    """
    from sortscore.utils.file_utils import ensure_output_subdirs
    
    # Ensure output directories exist
    ensure_output_subdirs(output_dir)
    
    # Save normalized scores
    scores_dir = os.path.join(output_dir, 'scores')
    scores_file = os.path.join(scores_dir, "batch_dna_scores.csv")
    
    dna_batch_data = results['normalized_scores'].copy()
    if not dna_batch_data.empty:
        dna_batch_data.to_csv(scores_file, index=False)
        logger.info(f"Saved batch DNA scores to {scores_file}")

    aa_batch_data = results.get('normalized_aa_scores')
    if aa_batch_data is None:
        aa_batch_data = pd.DataFrame()
    if not aa_batch_data.empty:
        aa_scores_file = os.path.join(scores_dir, "batch_aa_scores.csv")
        aa_batch_data.to_csv(aa_scores_file, index=False)
        logger.info(f"Saved batch AA scores to {aa_scores_file}")
    
    # Save combined statistics
    stats_file = os.path.join(scores_dir, "batch_stats.json")
    with open(stats_file, 'w') as f:
        json.dump(results['combined_stats'], f, indent=2)
    logger.info(f"Saved batch statistics to {stats_file}")
