"""
Heatmap visualization functions for Sort-seq variant analysis.

This module provides specialized heatmap plotting functions for MAVE-style visualizations,
including biophysical properties panels and codon-level analysis features.

Examples
--------
>>> from sortscore.visualization.heatmaps import plot_heatmap
>>> plot_heatmap(data, 'avgscore', experiment_config)
"""

import logging
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, List
from sortscore.visualization.heatmap_matrix import make_dms_matrix, fill_wt, make_col_avg_df, get_dropout
from sortscore.visualization.plots import generate_position_avg_colors
from sortscore.utils.load_experiment import ExperimentConfig
from sortscore.analysis.batch_config import BatchConfig
from sortscore.utils.sequence_parsing import convert_aa_to_three_letter


def _add_biophysical_properties_panel(ax, row_labels, aa_boundaries, is_small_heatmap=False):
    """Add biophysical properties heatmap panel for amino acid groups in codon heatmaps."""
    
    # Amino acid properties with numeric encoding for heatmap (both single and three-letter codes)
    aa_properties = {
        'W': {'type': 4}, 'Trp': {'type': 4},  # Aromatic
        'F': {'type': 4}, 'Phe': {'type': 4},  # Aromatic  
        'Y': {'type': 4}, 'Tyr': {'type': 4},  # Aromatic
        'P': {'type': 5}, 'Pro': {'type': 5},  # Special
        'M': {'type': 3}, 'Met': {'type': 3},  # Sulfur
        'I': {'type': 2}, 'Ile': {'type': 2},  # Branched
        'L': {'type': 2}, 'Leu': {'type': 2},  # Branched
        'V': {'type': 2}, 'Val': {'type': 2},  # Branched
        'A': {'type': 0}, 'Ala': {'type': 0},  # Small
        'G': {'type': 0}, 'Gly': {'type': 0},  # Small
        'C': {'type': 3}, 'Cys': {'type': 3},  # Sulfur
        'S': {'type': 1}, 'Ser': {'type': 1},  # Polar
        'T': {'type': 1}, 'Thr': {'type': 1},  # Polar
        'Q': {'type': 1}, 'Gln': {'type': 1},  # Polar
        'N': {'type': 1}, 'Asn': {'type': 1},  # Polar
        'D': {'type': 6}, 'Asp': {'type': 6},  # Acidic
        'E': {'type': 6}, 'Glu': {'type': 6},  # Acidic
        'H': {'type': 7}, 'His': {'type': 7},  # Basic
        'R': {'type': 7}, 'Arg': {'type': 7},  # Basic
        'K': {'type': 7}, 'Lys': {'type': 7},  # Basic
        '*': {'type': 8}, 'Ter': {'type': 8}   # Stop
    }

    # Pastel color mappings (only for type now) - Start uses same color as Sulfur
    type_colors = ['#F5F5F5', '#C8E6C8', '#FFD4B3', '#FFF2CC', 
                   '#E6D4FF', '#FFCCCC', '#CCE0FF', '#F4CCCC', '#D3D3D3', '#FFF2CC']
    type_labels = {0: 'Small', 1: 'Polar', 2: 'Branched', 3: 'Sulfur', 
                   4: 'Aromatic', 5: 'Special', 6: 'Acidic (–)', 7: 'Basic (+)', 8: 'Stop', 9: 'Start'}
    
    # Create property array for all rows
    type_array = np.zeros((len(row_labels), 1))
    
    for i, label in enumerate(row_labels):
        # Extract amino acid from label like "W(TGG)" -> "W" or just "W"
        aa = label.split('(')[0] if '(' in label else label
        if aa in aa_properties:
            type_array[i, 0] = aa_properties[aa]['type']
    
    # Plot only type column (charge information is now in the labels)
    for i in range(len(row_labels)):
        type_val = int(type_array[i, 0])
        
        # Add type rectangle - double width for small heatmaps, regular width for codon heatmaps
        width = 2 if is_small_heatmap else 1
        ax.add_patch(plt.Rectangle((0, i), width, 1, facecolor=type_colors[type_val],
                                  edgecolor=type_colors[type_val], linewidth=0))
    
    # Add labels for type blocks
    # Find continuous blocks of same type
    type_blocks = []
    current_type = type_array[0, 0]
    start_idx = 0
    for i in range(1, len(row_labels)):
        if type_array[i, 0] != current_type:
            type_blocks.append((int(current_type), start_idx, i-1))
            current_type = type_array[i, 0]
            start_idx = i
    # Add final block
    type_blocks.append((int(current_type), start_idx, len(row_labels)-1))
    
    # Add text labels for type blocks
    for type_val, start, end in type_blocks:
        center_y = (start + end + 1) / 2
        # Use horizontal rotation for both small and codon heatmaps
        rotation = 0
        # Use font size 16 for both small and codon heatmaps
        fontsize = 16
        # Position text in center of the (potentially wider) column
        x_center = 1 if is_small_heatmap else 0.5
        ax.text(x_center, center_y, type_labels[type_val], ha='center', va='center',
               fontsize=fontsize, fontweight='bold', color='black', rotation=rotation)
    
    # Set up the axis - adjust width based on heatmap type
    axis_width = 2 if is_small_heatmap else 1
    ax.set_xlim(0, axis_width)
    ax.set_ylim(0, len(row_labels))
    ax.invert_yaxis()
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


def plot_heatmap(
    data: pd.DataFrame,
    score_col: str,
    experiment: ExperimentConfig,
    wt_score: float = 1.0,
    fig_size: str = 'small',
    export_heatmap: bool = True,
    output: Optional[str] = None,
    fig_format: str = 'png',
    dpi: int = 300,
    tick_values: Optional[List[float]] = None,
    tick_labels: Optional[List[str]] = None,
    motif_indices: Optional[List[int]] = None,
    row_avg: bool = False,
    title: Optional[str] = None,
    export_matrix: bool = False,
    show_biophysical_properties: bool = False,
    export_positional_averages: bool = False,
    suffix: Optional[str] = None,
    heatmap_basename: str = "heatmap",
    matrix_basename: str = "heatmap_matrix",
    transparent: bool = True
) -> None:
    """
    Plot a MAVE heatmap using a matrix of activity scores.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    experiment : ExperimentConfig
        Experiment configuration dataclass instance.
    wt_score : float, default 1.0
        Score to assign to WT positions.
    fig_size : str, default 'small'
        Figure size ('small', 'large', 'long').
    export_heatmap : bool, default True
        If True, save the plot to file.
    output : str, optional
        Output directory for all exported files.
    fig_format : str, default 'png'
        Figure output format ('png' or 'svg').
    dpi : int, default 300
        Resolution for saved plot.
    tick_values : list of float, optional
        Tick values for the colorbar.
    tick_labels : list of str, optional
        Tick labels for the colorbar.
    motif_indices : list of int, optional
        Indices to highlight with transparency.
    row_avg : bool, default False
        If True, plot row averages.
    title : str, optional
        Plot title.
    export_matrix : bool, default False
        If True, export the heatmap matrix data to CSV.
    show_biophysical_properties : bool, default False
        If True, show biophysical properties panel beside the heatmap.
    export_positional_averages : bool, default False
        If True, export positional averages with hex colors to CSV for protein structure visualization.
    suffix : str, optional
        Suffix to append to output filenames for consistent naming.
    heatmap_basename : str, optional
        Base name for exported heatmap image (default: ``"heatmap"``).
    matrix_basename : str, optional
        Base name for exported heatmap matrix CSV (default: ``"heatmap_matrix"``).
    transparent : bool, default True
        If True, save plot with transparent background. If False, use white background.
    """
    logger = logging.getLogger(__name__)
    
    # Build output prefix once for all exports
    if output:
        output_prefix = f"{output}/{experiment.experiment_name}"
    else:
        output_prefix = experiment.experiment_name
    
    # Create DMS matrix and prepare data
    dms_matrix = make_dms_matrix(
        data,
        score_col,
        experiment.num_positions,
        experiment.wt_seq,
        experiment.variant_type,
        experiment.mutagenesis_variants,
        experiment.position_type
    )
    dropout_num, dropout_percent = get_dropout(dms_matrix, experiment.mutagenesis_variants)
    heatmap_df = fill_wt(dms_matrix, wt_score)
    col_avg_df = make_col_avg_df(heatmap_df)
    
    # Set up colormap and normalization (needed for exports)
    min_val = heatmap_df.min().min()
    max_val = heatmap_df.max().max()
    norm = plt.Normalize(vmin=min_val, vmax=max_val)
    cmap = plt.cm.magma
    
    # Export heatmap score matrix
    if export_matrix:
        matrix_for_export = heatmap_df.copy()
        # Update column names to reflect true residue positions
        new_columns = [str(i) for i in range(experiment.min_pos, experiment.min_pos + experiment.num_aa)]
        matrix_for_export.columns = new_columns

        if suffix:
            matrix_output = f"{output_prefix}_{matrix_basename}_{suffix}.csv"
        else:
            matrix_output = f"{output_prefix}_{matrix_basename}.csv"
        matrix_for_export.to_csv(matrix_output)
        logger.info(f"Heatmap matrix saved to {matrix_output}")
    
    # Export positional averages with colors for protein structure visualizations
    if export_positional_averages:
        # Use the same normalization and colormap as the main heatmap for exact color matching
        averages_colors = generate_position_avg_colors(heatmap_df, colormap='magma', norm=norm, cmap=cmap)
        if suffix:
            averages_output = f"{output_prefix}_positional_averages_{suffix}.csv"
        else:
            averages_output = f"{output_prefix}_positional_averages.csv"
        averages_colors.to_csv(averages_output, index=False)
        logger.info(f"Positional averages with colors saved to {averages_output}")

    # Set up masks and figure parameters
    nan_mask = dms_matrix.isnull()
    wt_mask = dms_matrix == 'WT'
    
    # Determine if this is a codon heatmap and whether to show biophysical properties
    num_rows = len(dms_matrix.index)
    is_codon_heatmap = num_rows > 21  # More than 21 rows = codon heatmap
    is_small_aa_heatmap = num_rows == 21  # Exactly 21 rows = standard AA heatmap
    show_props = show_biophysical_properties and (is_codon_heatmap or is_small_aa_heatmap)
    
    # Configure figure size
    if fig_size == 'small':
        width = max(16.5, experiment.num_aa * 0.15)
        height = 30 if is_codon_heatmap else 12
        facecolor = 'none' if transparent else 'white'
        fig = plt.figure(figsize=(width, height), facecolor=facecolor)
        tick_freq = max(1, experiment.num_aa // 15)
    elif fig_size == 'large':
        width = max(30, experiment.num_aa * 0.25)
        height = 35 if is_codon_heatmap else 10
        facecolor = 'none' if transparent else 'white'
        fig = plt.figure(figsize=(width, height), facecolor=facecolor)
        tick_freq = max(1, experiment.num_aa // 25)
    elif fig_size == 'long':
        width = max(30, experiment.num_aa * 0.3)
        height = 45 if is_codon_heatmap else 25
        facecolor = 'none' if transparent else 'white'
        fig = plt.figure(figsize=(width, height), facecolor=facecolor)
        tick_freq = max(1, experiment.num_aa // 40)

    # Set up subplot layout based on row_avg and biophysical properties
    if row_avg:
        row_avg_df = pd.DataFrame(heatmap_df.mean(axis=1), columns=['Avg'])
        if is_codon_heatmap:
            avg_height_ratio = 0.75 
            main_height_ratio = 60
        else:
            avg_height_ratio = 1
            main_height_ratio = 45
            
        if show_props:
            # Add space for biophysical properties panel
            gs = GridSpec(2,5, width_ratios=[1, 35, 5, 1, 1], height_ratios=[avg_height_ratio, main_height_ratio], hspace=0.03, wspace=0.05)
            ax1 = fig.add_subplot(gs[0, 1])
            ax2 = fig.add_subplot(gs[1, 1])
            ax_props = fig.add_subplot(gs[1, 2])  # Biophysical properties panel
            cax = fig.add_subplot(gs[1, 3])
            ax3 = fig.add_subplot(gs[1, 0])
        else:
            gs = GridSpec(2,3, width_ratios=[1, 35, 1], height_ratios=[avg_height_ratio, main_height_ratio], hspace=0.03, wspace=0.03)
            ax1 = fig.add_subplot(gs[0, 1])
            ax2 = fig.add_subplot(gs[1, 1])
            cax = fig.add_subplot(gs[1, 2])
            ax3 = fig.add_subplot(gs[1, 0])
            ax_props = None
    else:
        if is_codon_heatmap:
            avg_height_ratio = 0.75
            main_height_ratio = 35
        else:
            avg_height_ratio = 1
            main_height_ratio = 20
            
        if show_props:
            # Add space for biophysical properties panel
            gs = GridSpec(2, 4, width_ratios=[35, 5, 1, 1], height_ratios=[avg_height_ratio, main_height_ratio], hspace=0.03, wspace=0.05)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax_props = fig.add_subplot(gs[1, 1])  # Biophysical properties panel
            cax = fig.add_subplot(gs[1, 2])
        else:
            gs = GridSpec(2, 2, width_ratios=[35, 1], height_ratios=[avg_height_ratio, main_height_ratio], hspace=0.03, wspace=0.03)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            cax = fig.add_subplot(gs[1, 1])
            ax_props = None

    # Create heatmaps
    sns.heatmap(col_avg_df, annot=False, cmap=cmap, cbar=False, ax=ax1, norm=norm)
    ax = sns.heatmap(heatmap_df, cmap=cmap, cbar=False, ax=ax2, norm=norm)
    
    if nan_mask.any().any():
        sns.heatmap(nan_mask, annot=False, cmap=cmap, cbar=False, mask=~nan_mask, ax=ax2)
    
    if row_avg:
        sns.heatmap(row_avg_df, cmap=cmap, cbar=False, ax=ax3, norm=norm)
        
    # Add codon-specific features and biophysical properties
    if is_codon_heatmap or show_props:
        row_labels = list(heatmap_df.index)
        aa_boundaries = []
        current_aa = None
        
        # Find amino acid boundaries (for codon heatmaps or for grouping in small AA heatmaps)
        for i, label in enumerate(row_labels):
            # Extract amino acid from label like "W(TGG)" -> "W" or just "W"
            aa = label.split('(')[0] if '(' in label else label
            if current_aa is not None and aa != current_aa:
                aa_boundaries.append(i)
            current_aa = aa
            
        # Add horizontal white lines at amino acid boundaries (only for codon heatmaps)
        if is_codon_heatmap:
            for boundary in aa_boundaries:
                # Position lines at amino acid boundaries (between different amino acid groups)
                ax2.axhline(y=boundary, color='white', linewidth=2, solid_capstyle='butt')
                if row_avg:
                    ax3.axhline(y=boundary, color='white', linewidth=2, solid_capstyle='butt')


        # Add biophysical properties panel when requested
        if ax_props is not None:
            _add_biophysical_properties_panel(ax_props, row_labels, aa_boundaries, is_small_aa_heatmap)

    # Handle NaN and WT visualization
    for i in range(len(nan_mask)):
        for j in range(len(nan_mask.columns)):
            if nan_mask.iloc[i, j]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, hatch='//', edgecolor='lightgray', facecolor='white'))
    
    wt_indices = np.where(wt_mask)
    ax2.scatter(wt_indices[1] + 0.5, wt_indices[0] + 0.5, color='white', s=50, alpha=1.0, 
                edgecolors='lightgray', linewidths=0.5)
    
    if motif_indices:
        pass  # Motif highlighting logic would go here

    # Set up colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    if tick_values is not None and tick_labels is not None:
        cbar = plt.colorbar(sm, cax=cax, ticks=tick_values)
        cbar.ax.set_yticklabels(tick_labels, fontsize=18)
        cbar.ax.set_yticks(tick_values)
    else:
        cbar = plt.colorbar(sm, cax=cax)

    # Configure axes labels and ticks
    if experiment.position_type == 'dna':
        x_labels = [str(i) for i in range(1, experiment.num_positions + 1)]
    else:
        x_labels = [str(i) for i in range(experiment.min_pos, experiment.min_pos + experiment.num_aa)]
    
    tick_indices = list(range(0, len(x_labels), tick_freq))
    if 0 not in tick_indices:
        tick_indices = [0] + tick_indices
    
    ax.set_xticks([i + 0.5 for i in tick_indices])
    ax.set_xticklabels([x_labels[i] for i in tick_indices], rotation=0)

    # Set plot title
    if dropout_num > 0:
        plot_title = title or f'MAVE Heatmap - Dropout {dropout_num} variants ({round(dropout_percent)}%)'
    else:
        plot_title = title or 'MAVE Heatmap'
        
    if fig_size == 'small':
        fig.suptitle(plot_title, fontsize=32, y=0.89)
    elif fig_size == 'large':
        fig.suptitle(plot_title, fontsize=34, y=0.91)
        
    fig.subplots_adjust(top=0.85, bottom=0.15)

    # Configure subplot labels and formatting
    ax1.set_ylabel('Avg', fontsize=20, labelpad=20, ha='center', rotation=0)
    label = ax1.yaxis.get_label()
    label.set_verticalalignment('center')
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # X-axis label depends on position_type
    if experiment.position_type == 'dna':
        ax2.set_xlabel('DNA Position', fontsize=28)
    else:
        ax2.set_xlabel('Residue Sequence Number', fontsize=28)
        
    ax2.set_ylabel('Variant Amino Acid', fontsize=28)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    plt.setp(ax2.get_yticklabels(), rotation=0, ha='right')
    
    if row_avg:
        ax3.set_ylabel('Variant Amino Acid', fontsize=28, labelpad=8, ha='center')
        ax3.set_yticks(ax2.get_yticks())
        ax3.set_yticklabels(ax2.get_yticklabels())
        ax3.tick_params(axis='y', which='major', labelsize=20)
        plt.setp(ax3.get_yticklabels(), rotation=0, ha='right')
        ax2.set_ylabel('')
        ax2.set_yticks([])
        ax3.set_xticks([])

    # Save or show plot
    if export_heatmap:
        if suffix:
            heatmap_output = f"{output_prefix}_{heatmap_basename}_{suffix}.{fig_format}"
        else:
            heatmap_output = f"{output_prefix}_{heatmap_basename}.{fig_format}"
        facecolor = 'none' if transparent else 'white'
        plt.savefig(heatmap_output, dpi=dpi, format=fig_format, facecolor=facecolor, edgecolor='none')
        logger.info(f"Heatmap plot saved to {heatmap_output}")
    else:
        plt.show()
        
    plt.close()


def plot_tiled_heatmap(
    batch_data: pd.DataFrame,
    score_col: str,
    batch_config: BatchConfig,
    experiments: List[ExperimentConfig],
    wt_score: Optional[float] = None,
    fig_size: str = 'large',
    export: bool = False,
    output: Optional[str] = None,
    format: str = 'png',
    dpi: int = 300,
    tick_values: Optional[List[float]] = None,
    tick_labels: Optional[List[str]] = None,
    title: Optional[str] = None,
    export_matrix: bool = False,
    show_biophysical_properties: bool = False,
    transparent: bool = True
) -> None:
    """
    Plot a tiled heatmap combining data from multiple experiments with proper position mapping.
    
    This function creates a unified heatmap visualization for batch-processed experiments,
    handling position gaps and proper coordinate mapping for tiled experimental designs.
    
    Parameters
    ----------
    batch_data : pd.DataFrame
        Combined normalized data from all experiments with 'batch' column
    score_col : str
        Name of the score column to visualize
    batch_config : BatchConfig
        Batch configuration containing global position parameters
    experiments : List[ExperimentConfig]
        List of individual experiment configurations for position mapping
    wt_score : Optional[float]
        Wild-type reference score for visualization reference
    fig_size : str, optional
        Size preset ('small', 'medium', 'large'), default 'large'
    export : bool, optional
        Whether to save the heatmap to file, default False
    output : Optional[str]
        Output file path if exporting
    format : str, optional
        Export format ('png', 'pdf', 'svg'), default 'png'
    dpi : int, optional
        Resolution for exported image, default 300
    tick_values : Optional[List[float]]
        Custom colorbar tick values
    tick_labels : Optional[List[str]]
        Custom colorbar tick labels
    title : Optional[str]
        Custom heatmap title
    export_matrix : bool, optional
        Whether to export the underlying matrix data, default False
    show_biophysical_properties : bool, optional
        Whether to show biophysical properties panel, default False
    transparent : bool, default True
        If True, save plot with transparent background. If False, use white background.
        
    Examples
    --------
    Plot combined heatmap from batch analysis results:
    
    >>> plot_tiled_heatmap(
    ...     batch_data=normalized_scores,
    ...     score_col='avgscore',
    ...     batch_config=config,
    ...     experiments=experiment_list,
    ...     export=True,
    ...     output='combined_heatmap.png'
    ... )
    """
    logging.info(f"Creating tiled heatmap for {len(experiments)} experiments")
    
    # Determine global position range
    if batch_config.global_min_pos is not None and batch_config.global_max_pos is not None:
        global_min_pos = batch_config.global_min_pos
        global_max_pos = batch_config.global_max_pos
    else:
        # Calculate from experiment configs
        global_min_pos = min(exp.min_pos for exp in experiments)
        global_max_pos = max(exp.max_pos for exp in experiments)
    
    # Create position mapping for each experiment
    position_mapping = {}
    experiment_ranges = {}
    
    for i, experiment in enumerate(experiments, 1):
        batch_name = f'experiment{i}'
        position_mapping[batch_name] = {
            'min_pos': experiment.min_pos,
            'max_pos': experiment.max_pos,
            'global_offset': experiment.min_pos - global_min_pos
        }
        experiment_ranges[batch_name] = list(range(experiment.min_pos, experiment.max_pos + 1))
    
    # Get mutagenesis variants from first experiment (assume consistent across experiments)
    mutagenesis_variants = experiments[0].mutagenesis_variants
    if mutagenesis_variants is None:
        mutagenesis_variants = ['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']
    
    # Create global matrix with NaN values
    global_positions = list(range(global_min_pos, global_max_pos + 1))
    global_matrix = pd.DataFrame(
        index=mutagenesis_variants,
        columns=global_positions,
        dtype=float
    )
    global_matrix.iloc[:, :] = np.nan
    
    # Fill matrix with data from each experiment
    for batch_name, batch_df in batch_data.groupby('batch'):
        if batch_name in position_mapping:
            exp_mapping = position_mapping[batch_name]
            exp_positions = experiment_ranges[batch_name]
            
            # Create matrix for this experiment
            exp_matrix = make_dms_matrix(
                batch_df, 
                score_col, 
                exp_positions, 
                mutagenesis_variants,
                fill_na=True
            )
            
            # Map experiment matrix to global matrix
            for pos in exp_positions:
                if pos in exp_matrix.columns and pos in global_matrix.columns:
                    global_matrix[pos] = exp_matrix[pos]
    
    # Handle position breaks if requested
    if batch_config.allow_position_breaks:
        # Find gaps in data coverage
        data_coverage = ~global_matrix.isna().all(axis=0)
        covered_positions = global_matrix.columns[data_coverage].tolist()
        
        if len(covered_positions) > 0:
            # Use only covered positions for visualization
            plot_matrix = global_matrix[covered_positions]
            position_labels = covered_positions
        else:
            plot_matrix = global_matrix
            position_labels = global_positions
    else:
        # Use full global range
        plot_matrix = global_matrix
        position_labels = global_positions
    
    # Set figure size based on matrix dimensions
    size_presets = {
        'small': (12, 8),
        'medium': (16, 10),
        'large': (20, 12)
    }
    figsize = size_presets.get(fig_size, size_presets['large'])
    
    # Adjust figure width based on number of positions
    n_positions = len(position_labels)
    if n_positions > 100:
        figsize = (max(20, n_positions * 0.2), figsize[1])
    
    # Create figure
    facecolor = 'none' if transparent else 'white'
    fig = plt.figure(figsize=figsize, facecolor=facecolor)
    
    if show_biophysical_properties:
        # Create gridspec with biophysical properties panel
        gs = GridSpec(1, 2, figure=fig, width_ratios=[0.98, 0.02], wspace=0.01)
        ax_heatmap = fig.add_subplot(gs[0, 0])
        ax_props = fig.add_subplot(gs[0, 1])
    else:
        ax_heatmap = fig.add_subplot(111)
    
    # Create heatmap
    im = ax_heatmap.imshow(
        plot_matrix.values,
        aspect='auto',
        cmap='RdBu_r',
        interpolation='nearest'
    )
    
    # Set ticks and labels
    ax_heatmap.set_xticks(range(len(position_labels)))
    ax_heatmap.set_xticklabels(position_labels, rotation=45, ha='right')
    ax_heatmap.set_yticks(range(len(mutagenesis_variants)))
    ax_heatmap.set_yticklabels(mutagenesis_variants)
    
    # Add experiment boundaries as vertical lines if multiple experiments
    if len(experiments) > 1 and batch_config.allow_position_breaks:
        for i, batch_name in enumerate(position_mapping.keys()):
            if i > 0:  # Don't add line before first experiment
                exp_mapping = position_mapping[batch_name]
                # Find the boundary position in the plot
                boundary_pos = exp_mapping['min_pos']
                if boundary_pos in position_labels:
                    boundary_idx = position_labels.index(boundary_pos) - 0.5
                    ax_heatmap.axvline(x=boundary_idx, color='white', linewidth=2)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heatmap, shrink=0.8)
    if tick_values and tick_labels:
        cbar.set_ticks(tick_values)
        cbar.set_ticklabels(tick_labels)
    
    # Add biophysical properties panel if requested
    if show_biophysical_properties:
        _add_biophysical_properties_panel(
            ax_props, 
            mutagenesis_variants, 
            None,  # No boundaries needed for simple panel
            is_small_heatmap=fig_size == 'small'
        )
    
    # Set title
    if title:
        ax_heatmap.set_title(title)
    else:
        method = batch_config.batch_normalization_method
        n_exp = len(experiments)
        ax_heatmap.set_title(f'Combined Activity Scores ({method} normalization, {n_exp} experiments)')
    
    # Set axis labels
    ax_heatmap.set_xlabel('Position')
    ax_heatmap.set_ylabel('Amino Acid')
    
    # Export matrix data if requested
    if export_matrix and output:
        matrix_output = output.replace(f'.{format}', '_matrix.csv')
        plot_matrix.to_csv(matrix_output)
        logging.info(f"Exported matrix data to {matrix_output}")
    
    # Export figure if requested
    if export and output:
        plt.tight_layout()
        facecolor = 'none' if transparent else 'white'
        plt.savefig(output, dpi=dpi, bbox_inches='tight', format=format, facecolor=facecolor, edgecolor='none')
        logging.info(f"Exported tiled heatmap to {output}")
        
    plt.close()
