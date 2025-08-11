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
from sortscore.analysis.load_experiment import ExperimentConfig


def _add_biophysical_properties_panel(ax, row_labels, aa_boundaries):
    """Add biophysical properties heatmap panel for amino acid groups in codon heatmaps."""
    
    # Amino acid properties with numeric encoding for heatmap (indices match color arrays)
    aa_properties = {
        'W': {'type': 4, 'charge': 1},  # Aromatic, Neutral
        'F': {'type': 4, 'charge': 1},  # Aromatic, Neutral  
        'Y': {'type': 4, 'charge': 1},  # Aromatic, Neutral
        'P': {'type': 5, 'charge': 1},  # Special, Neutral
        'M': {'type': 3, 'charge': 1},  # Sulfur, Neutral
        'I': {'type': 2, 'charge': 1},  # Branched, Neutral
        'L': {'type': 2, 'charge': 1},  # Branched, Neutral
        'V': {'type': 2, 'charge': 1},  # Branched, Neutral
        'A': {'type': 0, 'charge': 1},  # Small, Neutral
        'G': {'type': 0, 'charge': 1},  # Small, Neutral
        'C': {'type': 3, 'charge': 1},  # Sulfur, Neutral
        'S': {'type': 1, 'charge': 1},  # Polar, Neutral
        'T': {'type': 1, 'charge': 1},  # Polar, Neutral
        'Q': {'type': 1, 'charge': 1},  # Polar, Neutral
        'N': {'type': 1, 'charge': 1},  # Polar, Neutral
        'D': {'type': 6, 'charge': 0},  # Acidic, Negative
        'E': {'type': 6, 'charge': 0},  # Acidic, Negative
        'H': {'type': 7, 'charge': 2},  # Basic, Positive
        'R': {'type': 7, 'charge': 2},  # Basic, Positive
        'K': {'type': 7, 'charge': 2},  # Basic, Positive
        '*': {'type': 8, 'charge': 1}   # Stop, Neutral
    }

    # Pastel color mappings
    charge_colors = ['#B3D9FF', '#F0F0F0', '#FFB3BA']  # Light blue, light gray, light red
    charge_labels = {0: '—', 1: '', 2: '+'}  # Em dash, blank, plus
    type_colors = ['#F5F5F5', '#C8E6C8', '#FFD4B3', '#FFF2CC', 
                   '#E6D4FF', '#FFCCCC', '#CCE0FF', '#F4CCCC', '#D3D3D3']
    type_labels = {0: 'Small', 1: 'Polar', 2: 'Branched', 3: 'Sulfur', 
                   4: 'Aromatic', 5: 'Special', 6: 'Acidic', 7: 'Basic', 8: 'Stop'}
    
    # Create property arrays for all rows
    charge_array = np.zeros((len(row_labels), 1))
    type_array = np.zeros((len(row_labels), 1))
    
    for i, label in enumerate(row_labels):
        # Extract amino acid from label like "W(TGG)" -> "W"
        aa = label.split('(')[0]
        if aa in aa_properties:
            charge_array[i, 0] = aa_properties[aa]['charge']
            type_array[i, 0] = aa_properties[aa]['type']
    
    # Plot charge and type columns
    for i in range(len(row_labels)):
        charge_val = int(charge_array[i, 0])
        type_val = int(type_array[i, 0])
        
        # Add charge rectangle
        ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=charge_colors[charge_val], 
                                  edgecolor=charge_colors[charge_val], linewidth=0))
        
        # Add type rectangle  
        ax.add_patch(plt.Rectangle((1, i), 1, 1, facecolor=type_colors[type_val],
                                  edgecolor=type_colors[type_val], linewidth=0))
    
    # Add single rotated labels for each property group spanning multiple rows
    # Find continuous blocks of same properties
    charge_blocks = []
    type_blocks = []
    
    # Process charge blocks
    current_charge = charge_array[0, 0]
    start_idx = 0
    for i in range(1, len(row_labels)):
        if charge_array[i, 0] != current_charge:
            if charge_labels[int(current_charge)]:  # Only add label if not empty (neutral)
                charge_blocks.append((int(current_charge), start_idx, i-1))
            current_charge = charge_array[i, 0]
            start_idx = i
    # Add final block
    if charge_labels[int(current_charge)]:
        charge_blocks.append((int(current_charge), start_idx, len(row_labels)-1))
    
    # Process type blocks
    current_type = type_array[0, 0]
    start_idx = 0
    for i in range(1, len(row_labels)):
        if type_array[i, 0] != current_type:
            type_blocks.append((int(current_type), start_idx, i-1))
            current_type = type_array[i, 0]
            start_idx = i
    # Add final block
    type_blocks.append((int(current_type), start_idx, len(row_labels)-1))
    
    # Add text labels for charge blocks
    for charge_val, start, end in charge_blocks:
        center_y = (start + end + 1) / 2
        ax.text(0.5, center_y, charge_labels[charge_val], ha='center', va='center',
               fontsize=18, fontweight='bold', color='black', rotation=0)
    
    # Add rotated text labels for type blocks
    for type_val, start, end in type_blocks:
        center_y = (start + end + 1) / 2
        ax.text(1.5, center_y, type_labels[type_val], ha='center', va='center',
               fontsize=18, fontweight='bold', color='black', rotation=90)
    
    # Set up the axis
    ax.set_xlim(0, 2)
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
    export: bool = False,
    output: Optional[str] = None,
    format: str = 'png',
    dpi: int = 300,
    tick_values: Optional[List[float]] = None,
    tick_labels: Optional[List[str]] = None,
    motif_indices: Optional[List[int]] = None,
    row_avg: bool = False,
    title: Optional[str] = None,
    export_matrix: bool = False
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
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    format : str, default 'png'
        Output format ('png' or 'svg').
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
    """
    logger = logging.getLogger(__name__)
    
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
    
    # Export matrix data if requested
    if export_matrix and export and output:
        # Create matrix with proper column headers (using experiment.min_pos)
        matrix_for_export = heatmap_df.copy()
        # Update column names to reflect true residue positions
        new_columns = [str(i) for i in range(experiment.min_pos, experiment.min_pos + experiment.num_aa)]
        matrix_for_export.columns = new_columns
        
        # Generate matrix output filename from plot output
        matrix_output = output.replace('.png', '_matrix.csv').replace('.svg', '_matrix.csv')
        matrix_for_export.to_csv(matrix_output)
        logger.info(f"Heatmap matrix saved to {matrix_output}")

    # Set up masks and figure parameters
    nan_mask = dms_matrix.isnull()
    wt_mask = dms_matrix == 'WT'
    
    # Determine if this is a codon heatmap
    num_rows = len(dms_matrix.index)
    is_codon_heatmap = num_rows > 21  # More than 21 rows = codon heatmap
    
    # Configure figure size
    if fig_size == 'small':
        width = max(16.5, experiment.num_aa * 0.15)
        height = 30 if is_codon_heatmap else 12
        fig = plt.figure(figsize=(width, height), facecolor='white')
        tick_freq = max(1, experiment.num_aa // 15)
    elif fig_size == 'large':
        width = max(30, experiment.num_aa * 0.25)
        height = 35 if is_codon_heatmap else 10
        fig = plt.figure(figsize=(width, height), facecolor='white')
        tick_freq = max(1, experiment.num_aa // 25)
    elif fig_size == 'long':
        width = max(30, experiment.num_aa * 0.3)
        height = 45 if is_codon_heatmap else 25
        fig = plt.figure(figsize=(width, height), facecolor='white')
        tick_freq = max(1, experiment.num_aa // 40)

    # Set up subplot layout based on row_avg and codon heatmap
    if row_avg:
        row_avg_df = pd.DataFrame(heatmap_df.mean(axis=1), columns=['Avg'])
        if is_codon_heatmap:
            avg_height_ratio = 0.75 
            main_height_ratio = 60
        else:
            avg_height_ratio = 1
            main_height_ratio = 45
            
        if is_codon_heatmap:
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
            
        if is_codon_heatmap:
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

    # Set up colormap and normalization
    min_val = heatmap_df.min().min()
    max_val = heatmap_df.max().max()
    norm = plt.Normalize(vmin=min_val, vmax=max_val)
    cmap = plt.cm.magma
    
    # Create heatmaps
    sns.heatmap(col_avg_df, annot=False, cmap=cmap, cbar=False, ax=ax1, norm=norm)
    ax = sns.heatmap(heatmap_df, cmap=cmap, cbar=False, ax=ax2, norm=norm)
    
    if nan_mask.any().any():
        sns.heatmap(nan_mask, annot=False, cmap=cmap, cbar=False, mask=~nan_mask, ax=ax2)
    
    if row_avg:
        sns.heatmap(row_avg_df, cmap=cmap, cbar=False, ax=ax3, norm=norm)
        
    # Add codon-specific features
    if is_codon_heatmap:
        row_labels = list(heatmap_df.index)
        aa_boundaries = []
        current_aa = None
        
        for i, label in enumerate(row_labels):
            # Extract amino acid from label like "W(TGG)" -> "W"
            aa = label.split('(')[0]
            if current_aa is not None and aa != current_aa:
                aa_boundaries.append(i)
            current_aa = aa
            
        # Add horizontal white lines at amino acid boundaries
        for boundary in aa_boundaries:
            ax2.axhline(y=boundary, color='white', linewidth=1)
            if row_avg:
                ax3.axhline(y=boundary, color='white', linewidth=1)
                
        # Add thick white dotted line below M(ATG) to separate M from C within sulfur group
        for i, label in enumerate(row_labels):
            aa = label.split('(')[0]
            if aa == 'M':  # Found M(ATG) row
                ax2.axhline(y=i+1, color='white', linewidth=3, linestyle=(0, (5, 2)))
                if row_avg:
                    ax3.axhline(y=i+1, color='white', linewidth=3, linestyle=(0, (5, 2)))
                
        # Add biophysical properties panel for codon heatmaps
        if ax_props is not None:
            _add_biophysical_properties_panel(ax_props, row_labels, aa_boundaries)

    # Handle NaN and WT visualization
    for i in range(len(nan_mask)):
        for j in range(len(nan_mask.columns)):
            if nan_mask.iloc[i, j]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, hatch='//', edgecolor='lightgray', facecolor='white'))
    
    wt_indices = np.where(wt_mask)
    ax2.scatter(wt_indices[1] + 0.5, wt_indices[0] + 0.5, color='white', s=30, alpha=0.5)
    
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
    if export and output:
        plt.savefig(output, dpi=dpi, format=format, facecolor='white', edgecolor='none')
        logger.info(f"Heatmap plot saved to {output}")
    else:
        plt.show()
        
    plt.close()