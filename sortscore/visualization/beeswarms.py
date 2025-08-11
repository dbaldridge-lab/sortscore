"""
Beeswarm and swarm plot visualization functions for Sort-seq variant analysis.

This module provides functions for creating beeswarm plots, violin plots with swarm overlays,
and grouped swarm plots for visualizing activity score distributions.

Examples
--------
>>> from sortscore.visualization.beeswarms import plot_beeswarm, plot_violin_swarm
>>> plot_beeswarm(df, x='annotation', y='avgscore')
"""

import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, List, Dict, Any
from matplotlib.patches import Rectangle

def plot_beeswarm(
    df: pd.DataFrame,
    x: str,
    y: str,
    hue: Optional[str] = None,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    size: float = 4,
    alpha: float = 0.7,
    dodge: bool = True
) -> None:
    """
    Create a beeswarm plot of activity scores.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis (e.g., annotation).
    y : str
        Column name for y-axis (score).
    hue : str, optional
        Column name for hue grouping.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    size : float, default 4
        Size of the data points.
    alpha : float, default 0.7
        Transparency of the data points.
    dodge : bool, default True
        If True, separate hue levels along categorical axis.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(15, 10))
    
    sns.swarmplot(data=df, x=x, y=y, hue=hue, size=size, alpha=alpha, 
                 dodge=dodge, legend=False)
    
    plt.title(title or 'Beeswarm Plot', fontsize=22)
    plt.xlabel(x, fontsize=18)
    plt.ylabel(y, fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    # Rotate x-axis labels if they're too long
    x_labels = df[x].unique()
    if any(len(str(label)) > 8 for label in x_labels):
        plt.xticks(rotation=45, ha='right')
    
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Beeswarm plot saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_violin_swarm(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    swarm_size: float = 3,
    swarm_alpha: float = 0.6
) -> None:
    """
    Create violin plot with overlaid swarm plot.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis grouping variable.
    y : str
        Column name for y-axis (score).
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    swarm_size : float, default 3
        Size of swarm plot points.
    swarm_alpha : float, default 0.6
        Transparency of swarm points.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(14, 10))
    
    # Create violin plot
    sns.violinplot(data=df, x=x, y=y, inner=None, alpha=0.8)
    
    # Overlay swarm plot
    sns.swarmplot(data=df, x=x, y=y, color='black', size=swarm_size, 
                 alpha=swarm_alpha)
    
    plt.title(title or 'Violin Plot with Swarm Overlay', fontsize=20)
    plt.xlabel(x, fontsize=16)
    plt.ylabel(y, fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Violin-swarm plot saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_strip_jitter(
    df: pd.DataFrame,
    x: str,
    y: str,
    hue: Optional[str] = None,
    jitter: float = 0.3,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    size: float = 5,
    alpha: float = 0.7
) -> None:
    """
    Create strip plot with jitter for categorical data visualization.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis grouping variable.
    y : str
        Column name for y-axis (score).
    hue : str, optional
        Column name for hue grouping.
    jitter : float, default 0.3
        Amount of jitter (spread) in x-direction.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    size : float, default 5
        Size of the data points.
    alpha : float, default 0.7
        Transparency of the data points.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 8))
    
    sns.stripplot(data=df, x=x, y=y, hue=hue, size=size, alpha=alpha,
                 jitter=jitter, dodge=True)
    
    plt.title(title or 'Strip Plot with Jitter', fontsize=18)
    plt.xlabel(x, fontsize=14)
    plt.ylabel(y, fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Strip plot saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_grouped_swarm(
    df: pd.DataFrame,
    x: str,
    y: str,
    hue: str,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    size: float = 4,
    alpha: float = 0.8,
    palette: Optional[str] = None
) -> None:
    """
    Create grouped swarm plot with multiple hue levels.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis grouping variable.
    y : str
        Column name for y-axis (score).
    hue : str
        Column name for hue grouping (creates separate groups).
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    size : float, default 4
        Size of the data points.
    alpha : float, default 0.8
        Transparency of the data points.
    palette : str, optional
        Color palette name for hue levels.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(16, 10))
    
    sns.swarmplot(data=df, x=x, y=y, hue=hue, size=size, alpha=alpha,
                 dodge=True, palette=palette)
    
    plt.title(title or f'Grouped Swarm Plot: {y} by {x} and {hue}', fontsize=20)
    plt.xlabel(x, fontsize=16)
    plt.ylabel(y, fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(title=hue, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', 
                   edgecolor='none', bbox_inches='tight')
        logger.info(f"Grouped swarm plot saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_boxplot_swarm(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    swarm_size: float = 3,
    swarm_alpha: float = 0.5,
    box_alpha: float = 0.7
) -> None:
    """
    Create box plot with overlaid swarm plot.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis grouping variable.
    y : str
        Column name for y-axis (score).
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    swarm_size : float, default 3
        Size of swarm plot points.
    swarm_alpha : float, default 0.5
        Transparency of swarm points.
    box_alpha : float, default 0.7
        Transparency of box plots.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 8))
    
    # Create box plot
    box_plot = sns.boxplot(data=df, x=x, y=y, showfliers=False, 
                          boxprops={'alpha': box_alpha})
    
    # Overlay swarm plot
    sns.swarmplot(data=df, x=x, y=y, color='black', size=swarm_size, 
                 alpha=swarm_alpha)
    
    plt.title(title or 'Box Plot with Swarm Overlay', fontsize=18)
    plt.xlabel(x, fontsize=14)
    plt.ylabel(y, fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Box-swarm plot saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_swarm_grid(
    df: pd.DataFrame,
    x: str,
    y: str,
    col: str,
    row: Optional[str] = None,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    size: float = 3
) -> None:
    """
    Create a grid of swarm plots using FacetGrid.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis within each facet.
    y : str
        Column name for y-axis (score).
    col : str
        Column name for faceting columns.
    row : str, optional
        Column name for faceting rows.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    size : float, default 3
        Size of the data points.
    """
    logger = logging.getLogger(__name__)
    
    # Create FacetGrid
    if row:
        g = sns.FacetGrid(df, col=col, row=row, height=5, aspect=1.2)
    else:
        g = sns.FacetGrid(df, col=col, height=5, aspect=1.2, col_wrap=3)
    
    # Map swarm plot to each facet
    g.map(sns.swarmplot, x, y, size=size, alpha=0.7)
    
    # Customize each subplot
    for ax in g.axes.flat:
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Set labels
    g.set_axis_labels(x, y)
    
    if title:
        g.fig.suptitle(title, fontsize=16, y=1.02)
    
    plt.tight_layout()
    
    if export and output:
        g.savefig(output, dpi=300, format='png', facecolor='white', 
                 edgecolor='none', bbox_inches='tight')
        logger.info(f"Swarm grid saved to {output}")
    else:
        plt.show()
    
    plt.close()


def plot_swarm_with_stats(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    show_mean: bool = True,
    show_median: bool = True,
    size: float = 4
) -> Dict[str, Dict[str, float]]:
    """
    Create swarm plot with overlaid statistical markers.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    x : str
        Column name for x-axis grouping variable.
    y : str
        Column name for y-axis (score).
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    show_mean : bool, default True
        If True, show mean markers.
    show_median : bool, default True
        If True, show median markers.
    size : float, default 4
        Size of swarm points.
        
    Returns
    -------
    dict
        Dictionary containing group statistics.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 8))
    
    # Create swarm plot
    sns.swarmplot(data=df, x=x, y=y, size=size, alpha=0.7)
    
    # Calculate and overlay statistics
    group_stats = {}
    x_positions = range(len(df[x].unique()))
    groups = sorted(df[x].unique())
    
    for i, group in enumerate(groups):
        group_data = df[df[x] == group][y].dropna()
        
        if len(group_data) == 0:
            continue
            
        mean_val = group_data.mean()
        median_val = group_data.median()
        
        group_stats[group] = {
            'mean': mean_val,
            'median': median_val,
            'std': group_data.std(),
            'n': len(group_data)
        }
        
        if show_mean:
            plt.scatter(i, mean_val, color='red', s=100, marker='D', 
                       label='Mean' if i == 0 else "", zorder=5)
        
        if show_median:
            plt.scatter(i, median_val, color='blue', s=100, marker='s',
                       label='Median' if i == 0 else "", zorder=5)
    
    plt.title(title or f'Swarm Plot with Statistics: {y} by {x}', fontsize=16)
    plt.xlabel(x, fontsize=14)
    plt.ylabel(y, fontsize=14)
    plt.xticks(rotation=45, ha='right')
    
    if show_mean or show_median:
        plt.legend()
    
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Swarm plot with stats saved to {output}")
    else:
        plt.show()
    
    plt.close()
    return group_stats