"""Utilities for generating heatmap visualizations from analysis outputs.

This module centralizes the heatmap generation logic. It prepares the data for amino acid
and codon heatmaps and delegates plotting to `plot_heatmap`.
"""

from __future__ import annotations

import logging
import os
from types import SimpleNamespace
from typing import Optional, Tuple, List

import pandas as pd

from sortscore.analysis.data_processing import aggregate_aa_data
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.visualization.heatmaps import plot_heatmap


def _score_column_name(experiment: ExperimentConfig) -> str:
    """Return the appropriate score column name based on avg_method."""
    if experiment.avg_method == 'simple-avg':
        return 'avgscore'
    suffix = experiment.avg_method.replace('-', '_')
    return f'avgscore_{suffix}'


def _build_aa_heatmap_config(experiment: ExperimentConfig) -> SimpleNamespace:
    """Create a lightweight config object tailored for AA heatmaps."""
    # AA heatmaps always use AA coordinates and an AA wild-type sequence.
    # The plotting code expects an "experiment-like" object with the same
    # attributes as ExperimentConfig (including num_positions).
    from sortscore.sequence_parsing import get_reference_sequence

    wt_aa_seq = get_reference_sequence(experiment.wt_seq, "aa")
    return SimpleNamespace(
        experiment_name=experiment.experiment_name,
        num_aa=experiment.num_aa,
        num_positions=experiment.num_aa,
        min_pos=experiment.min_pos,
        max_pos=experiment.max_pos,
        wt_seq=wt_aa_seq,
        variant_type='aa',
        mutagenesis_variants=experiment.mutagenesis_variants,
        position_type='aa',
    )


def _compute_wt_score(
    data: pd.DataFrame,
    score_col: str,
    annotate_col: str,
    annotate_value: str
) -> pd.NA | float:
    """Compute WT score from annotated rows if present."""
    if annotate_col not in data.columns:
        return pd.NA
    subset = data[data[annotate_col] == annotate_value]
    if len(subset) == 0 or score_col not in subset.columns:
        return pd.NA
    score_val = subset[score_col].mean()
    return float(score_val) if pd.notna(score_val) else pd.NA


def _tick_marks(
    values: pd.Series,
    wt_score: Optional[float],
    wt_label: str
) -> Tuple[Optional[List[float]], Optional[List[str]]]:
    """Generate tick positions and labels for colorbars."""
    if wt_score is None or pd.isna(wt_score):
        return None, None
    data_min = values.min()
    data_max = values.max()
    if pd.isna(data_min) or pd.isna(data_max):
        return None, None
    tick_values = [data_min, wt_score, data_max]
    tick_labels = [f'{data_min:.0f}', wt_label, f'{data_max:.0f}']
    return tick_values, tick_labels


def generate_heatmap_visualizations(
    scores_df: pd.DataFrame,
    experiment: ExperimentConfig,
    output_dir: str,
    output_suffix: str,
    fig_format: str = 'png',
    export_positional_averages: bool = False,
) -> None:
    """
    Generate amino acid and codon heatmaps for a completed analysis run.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing final scored variants.
    experiment : ExperimentConfig
        Experiment configuration used for plotting context.
    output_dir : str
        Base output directory (figures will be written under ``output_dir/figures``).
    output_suffix : str
        Suffix applied to exported file names for traceability.
    fig_format : str, optional
        Export figure format (default: ``'png'``).
    export_positional_averages : bool, optional
        Whether to export positional averages with colors for structure visualization.
    """
    logger = logging.getLogger(__name__)
    logger.info("Generating visualizations...")

    score_col = _score_column_name(experiment)
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Amino acid heatmap generation
    if 'aa_seq_diff' in scores_df.columns:
        if experiment.variant_type == 'aa':
            logger.info("AA data creation: scores_df shape %s", scores_df.shape)
            aa_data = scores_df[['aa_seq_diff', 'annotate_aa', score_col]].copy()
            logger.info("AA data created: aa_data shape %s", aa_data.shape)
        else:
            # For DNA variant type, aggregate to AA level
            aa_data = aggregate_aa_data(scores_df, score_col)
            logger.info(
                "AA data aggregated: aa_data shape %s; columns: %s",
                aa_data.shape,
                ", ".join(aa_data.columns),
            )

        wt_score = _compute_wt_score(aa_data, score_col, 'annotate_aa', 'synonymous')
        if pd.notna(wt_score):
            logger.info("Synonymous-WT score from AA data: %s", wt_score)
        else:
            logger.info("No synonymous WT score found for AA heatmap")

        tick_values, tick_labels = _tick_marks(aa_data[score_col], float(wt_score) if pd.notna(wt_score) else None, "WT Avg")
        aa_config = _build_aa_heatmap_config(experiment)

        plot_kwargs = dict(
            tick_values=tick_values,
            tick_labels=tick_labels,
            export_heatmap=True,
            output=figures_dir,
            fig_format=fig_format,
            export_matrix=True,
            export_positional_averages=export_positional_averages,
            show_biophysical_properties=experiment.biophysical_prop,
            suffix=output_suffix,
        )

        if pd.notna(wt_score):
            plot_heatmap(aa_data, score_col, aa_config, wt_score=float(wt_score), **plot_kwargs)
        else:
            plot_heatmap(aa_data, score_col, aa_config, **plot_kwargs)

        aa_heatmap_file = os.path.join(figures_dir, f"aa_heatmap_{experiment.avg_method}_{output_suffix}.{fig_format}")
        logger.info("Saved AA heatmap to %s", aa_heatmap_file)

    # Codon heatmap generation
    if experiment.variant_type == 'dna':
        wt_score_codon = _compute_wt_score(scores_df, score_col, 'annotate_dna', 'wt_dna')
        if pd.notna(wt_score_codon):
            logger.info("Found WT score for codon heatmap: %s", wt_score_codon)
        else:
            logger.info("No WT score found for codon heatmap")

        tick_values_codon, tick_labels_codon = _tick_marks(
            scores_df[score_col],
            float(wt_score_codon) if pd.notna(wt_score_codon) else None,
            "WT",
        )

        plot_kwargs = dict(
            tick_values=tick_values_codon,
            tick_labels=tick_labels_codon,
            export_heatmap=True,
            output=figures_dir,
            fig_format=fig_format,
            export_matrix=True,
            export_positional_averages=export_positional_averages,
            show_biophysical_properties=experiment.biophysical_prop,
            suffix=output_suffix,
        )

        if pd.notna(wt_score_codon):
            plot_heatmap(scores_df, score_col, experiment, wt_score=float(wt_score_codon), **plot_kwargs)
        else:
            plot_heatmap(scores_df, score_col, experiment, **plot_kwargs)

        codon_heatmap_file = os.path.join(figures_dir, f"codon_heatmap_{experiment.avg_method}_{output_suffix}.{fig_format}")
        logger.info("Saved codon heatmap to %s", codon_heatmap_file)
