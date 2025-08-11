"""
Histogram visualization functions for Sort-seq variant analysis.

This module provides functions for creating various types of histograms to analyze
activity score distributions, including grouped histograms, stacked histograms, and
comparative histogram analyses.

Examples
--------
>>> from sortscore.visualization.histograms import plot_histogram, plot_stacked_histogram
>>> plot_histogram(df, 'avgscore', group_col='annotation')
"""

import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, List, Dict, Any
from matplotlib.patches import Rectangle

def plot_histogram(
    df: pd.DataFrame,
    score_col: str,
    group_col: str = 'annotate_dna',
    export: bool = False,
    output: Optional[str] = None,
    bins: int = 50,
    alpha: float = 0.7,
    highlight_groups: Optional[List[str]] = None
) -> None:
    """
    Plot histograms of activity scores grouped by annotation.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    group_col : str, default 'annotate_dna'
        Column name for grouping.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    bins : int, default 50
        Number of bins for histograms.
    alpha : float, default 0.7
        Transparency of histogram bars.
    highlight_groups : list of str, optional
        Groups to highlight with different styling.
    """
    logger = logging.getLogger(__name__)
    
    min_score = df[score_col].min()
    max_score = df[score_col].max()
    bin_edges = np.linspace(min_score, max_score, bins + 1)
    
    unique_values = df[group_col].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_values)))
    
    # Default highlighted groups
    if highlight_groups is None:
        highlight_groups = ['spike-in', 'wt_dna']
    
    other_values = [v for v in unique_values if v not in highlight_groups]
    
    plt.figure(figsize=(12, 8))
    
    # Plot non-highlighted groups first
    for value in other_values:
        color = colors[list(unique_values).index(value)]
        subset = df[df[group_col] == value]
        plt.hist(subset[score_col], bins=bin_edges, color=color, alpha=alpha, 
                label=f"{str(value)} (n={len(subset)})")
    
    # Plot highlighted groups with special styling
    for value in highlight_groups:
        if value in unique_values:
            color = colors[list(unique_values).index(value)]
            subset = df[df[group_col] == value]
            plt.hist(subset[score_col], bins=bin_edges, color=color, alpha=0.9, 
                    label=f"{str(value)} (n={len(subset)})", edgecolor='black', linewidth=1.5)
    
    plt.title(f'Histogram of {score_col} by {group_col}', fontsize=16)
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(title=group_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=1200, format='png', transparent=True, bbox_inches='tight')
        logger.info(f"Histogram plot saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_stacked_histogram(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    bins: int = 30,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    normalize: bool = False
) -> None:
    """
    Create stacked histogram showing group contributions to each bin.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    group_col : str
        Column name for grouping variable.
    bins : int, default 30
        Number of bins for histogram.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    normalize : bool, default False
        If True, normalize each bin to show proportions.
    """
    logger = logging.getLogger(__name__)
    
    # Prepare data for stacking
    groups = sorted(df[group_col].unique())
    colors = plt.cm.Set3(np.linspace(0, 1, len(groups)))
    
    # Calculate bin edges
    min_score = df[score_col].min()
    max_score = df[score_col].max()
    bin_edges = np.linspace(min_score, max_score, bins + 1)
    
    plt.figure(figsize=(12, 8))
    
    # Collect data for each group
    group_data = []
    group_labels = []
    
    for group in groups:
        group_scores = df[df[group_col] == group][score_col].dropna()
        if len(group_scores) > 0:
            group_data.append(group_scores)
            group_labels.append(f"{group} (n={len(group_scores)})")
    
    # Create stacked histogram
    plt.hist(group_data, bins=bin_edges, label=group_labels, stacked=True,
            alpha=0.8, color=colors[:len(group_data)], density=normalize)
    
    plt.title(title or f'Stacked Histogram: {score_col} by {group_col}', fontsize=16)
    plt.xlabel('Activity Score', fontsize=14)
    
    if normalize:
        plt.ylabel('Density', fontsize=14)
    else:
        plt.ylabel('Frequency', fontsize=14)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', 
                   edgecolor='none', bbox_inches='tight')
        logger.info(f"Stacked histogram saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_histogram_grid(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    bins: int = 20,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    ncols: int = 3,
    figsize_per_subplot: tuple = (4, 3)
) -> None:
    """
    Create a grid of histograms, one for each group.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    group_col : str
        Column name for grouping variable.
    bins : int, default 20
        Number of bins per histogram.
    title : str, optional
        Overall plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    ncols : int, default 3
        Number of columns in the grid.
    figsize_per_subplot : tuple, default (4, 3)
        Size of each subplot.
    """
    logger = logging.getLogger(__name__)
    
    groups = sorted(df[group_col].unique())
    ngroups = len(groups)
    nrows = int(np.ceil(ngroups / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, 
                           figsize=(figsize_per_subplot[0] * ncols, 
                                   figsize_per_subplot[1] * nrows))
    
    if ngroups == 1:
        axes = [axes]
    elif nrows == 1:
        axes = axes if isinstance(axes, list) else [axes]
    else:
        axes = axes.flatten()
    
    # Common bin range for all subplots
    min_score = df[score_col].min()
    max_score = df[score_col].max()
    bin_edges = np.linspace(min_score, max_score, bins + 1)
    
    colors = plt.cm.tab10(np.linspace(0, 1, ngroups))
    
    for i, group in enumerate(groups):
        ax = axes[i]
        group_data = df[df[group_col] == group][score_col].dropna()
        
        ax.hist(group_data, bins=bin_edges, color=colors[i], alpha=0.7, 
               edgecolor='black', linewidth=0.5)
        ax.set_title(f"{group} (n={len(group_data)})", fontsize=12)
        ax.set_xlabel('Activity Score', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Add mean line
        if len(group_data) > 0:
            mean_val = group_data.mean()
            ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, alpha=0.8)
    
    # Hide unused subplots
    for i in range(ngroups, len(axes)):
        axes[i].set_visible(False)
    
    if title:
        fig.suptitle(title, fontsize=16)
    
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', 
                   edgecolor='none', bbox_inches='tight')
        logger.info(f"Histogram grid saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_normalized_comparison(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    reference_group: Optional[str] = None,
    bins: int = 30,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None
) -> None:
    """
    Plot normalized histograms for direct comparison of distributions.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    group_col : str
        Column name for grouping variable.
    reference_group : str, optional
        Group to highlight as reference.
    bins : int, default 30
        Number of bins for histograms.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    """
    logger = logging.getLogger(__name__)
    
    groups = sorted(df[group_col].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
    
    plt.figure(figsize=(12, 8))
    
    # Common bin range
    min_score = df[score_col].min()
    max_score = df[score_col].max()
    bin_edges = np.linspace(min_score, max_score, bins + 1)
    
    for i, group in enumerate(groups):
        group_data = df[df[group_col] == group][score_col].dropna()
        
        if len(group_data) == 0:
            continue
        
        # Special styling for reference group
        if group == reference_group:
            alpha = 0.9
            linewidth = 2
            linestyle = '-'
            color = 'black'
        else:
            alpha = 0.7
            linewidth = 1
            linestyle = '-'
            color = colors[i]
        
        plt.hist(group_data, bins=bin_edges, density=True, alpha=alpha, 
                color=color, label=f"{group} (n={len(group_data)})",
                histtype='step', linewidth=linewidth, linestyle=linestyle)
    
    plt.title(title or f'Normalized Distribution Comparison: {score_col}', fontsize=16)
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', 
                   edgecolor='none', bbox_inches='tight')
        logger.info(f"Normalized comparison saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_cumulative_histogram(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    bins: int = 100
) -> None:
    """
    Plot cumulative distribution functions for different groups.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    group_col : str
        Column name for grouping variable.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    bins : int, default 100
        Number of bins for CDF calculation.
    """
    logger = logging.getLogger(__name__)
    
    groups = sorted(df[group_col].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
    
    plt.figure(figsize=(12, 8))
    
    for i, group in enumerate(groups):
        group_data = df[df[group_col] == group][score_col].dropna()
        
        if len(group_data) == 0:
            continue
        
        # Calculate cumulative distribution
        sorted_data = np.sort(group_data)
        cumulative = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        
        plt.plot(sorted_data, cumulative, color=colors[i], linewidth=2,
                label=f"{group} (n={len(group_data)})")
    
    plt.title(title or f'Cumulative Distribution: {score_col}', fontsize=16)
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Cumulative Probability', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.xlim(df[score_col].min(), df[score_col].max())
    plt.ylim(0, 1)
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', 
                   edgecolor='none', bbox_inches='tight')
        logger.info(f"Cumulative histogram saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_histogram_with_stats(
    df: pd.DataFrame,
    score_col: str,
    bins: int = 50,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    show_normal_overlay: bool = True,
    show_quartiles: bool = True
) -> Dict[str, float]:
    """
    Plot histogram with statistical overlays and annotations.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    bins : int, default 50
        Number of bins for histogram.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    show_normal_overlay : bool, default True
        If True, overlay fitted normal distribution.
    show_quartiles : bool, default True
        If True, show quartile lines.
        
    Returns
    -------
    dict
        Dictionary of calculated statistics.
    """
    logger = logging.getLogger(__name__)
    
    clean_data = df[score_col].dropna()
    
    # Calculate statistics
    stats_dict = {
        'mean': clean_data.mean(),
        'median': clean_data.median(),
        'std': clean_data.std(),
        'q25': clean_data.quantile(0.25),
        'q75': clean_data.quantile(0.75),
        'n': len(clean_data)
    }
    
    plt.figure(figsize=(12, 8))
    
    # Create histogram
    n, bins_used, patches = plt.hist(clean_data, bins=bins, density=True, 
                                    alpha=0.7, color='skyblue', edgecolor='black')
    
    # Add normal overlay
    if show_normal_overlay:
        x = np.linspace(clean_data.min(), clean_data.max(), 200)
        normal_pdf = (1 / (stats_dict['std'] * np.sqrt(2 * np.pi))) * \
                    np.exp(-0.5 * ((x - stats_dict['mean']) / stats_dict['std'])**2)
        plt.plot(x, normal_pdf, 'r-', linewidth=2, label='Normal fit')
    
    # Add statistical lines
    plt.axvline(stats_dict['mean'], color='red', linestyle='--', linewidth=2, 
               label=f"Mean: {stats_dict['mean']:.3f}")
    plt.axvline(stats_dict['median'], color='green', linestyle='--', linewidth=2,
               label=f"Median: {stats_dict['median']:.3f}")
    
    if show_quartiles:
        plt.axvline(stats_dict['q25'], color='orange', linestyle=':', linewidth=2,
                   label=f"Q1: {stats_dict['q25']:.3f}")
        plt.axvline(stats_dict['q75'], color='orange', linestyle=':', linewidth=2,
                   label=f"Q3: {stats_dict['q75']:.3f}")
    
    plt.title(title or f'Histogram with Statistics: {score_col}', fontsize=16)
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add statistics text box
    stats_text = (f"n = {stats_dict['n']}\n"
                 f"μ = {stats_dict['mean']:.3f}\n" 
                 f"σ = {stats_dict['std']:.3f}")
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
            verticalalignment='top', fontsize=12,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', 
                   edgecolor='none', bbox_inches='tight')
        logger.info(f"Histogram with stats saved to {output}")
    else:
        plt.show()
    
    plt.close()
    return stats_dict