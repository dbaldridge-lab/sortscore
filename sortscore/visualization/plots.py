"""
Visualization functions for Sort-seq variant analysis.

This module provides plotting functions for variant activity score data, including beeswarm and histogram plots.

Examples
--------
>>> from sortscore.visualization.plots import plot_activity_score_distribution, plot_beeswarm, plot_histogram
>>> from sortscore.visualization.heatmaps import plot_heatmap
"""
import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, List, Any
# Import plot_heatmap from the dedicated heatmaps module
from sortscore.visualization.heatmaps import plot_heatmap


def plot_activity_score_distribution(
    df: pd.DataFrame,
    score_col: str = 'avgscore',
    bins: int = 50,
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> None:
    """
    Plot the distribution of activity scores.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing activity scores.
    score_col : str, default 'avgscore'
        Column name for the activity score.
    bins : int, default 50
        Number of bins for the histogram.
    title : str, optional
        Title for the plot.
    save_path : str, optional
        If provided, save the plot to this file path.

    Returns
    -------
    None

    Examples
    --------
    >>> plot_activity_score_distribution(df, score_col='avgscore_rep_weighted')
    """
    logger = logging.getLogger(__name__)
    if score_col not in df.columns:
        logger.error(f"Score column '{score_col}' not found in DataFrame.")
        raise ValueError(f"Score column '{score_col}' not found in DataFrame.")
    plt.figure(figsize=(8, 5))
    plt.hist(df[score_col], bins=bins, color='skyblue', edgecolor='black')
    plt.xlabel('Activity Score')
    plt.ylabel('Count')
    plt.title(title or f'Activity Score Distribution ({score_col})')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
        logger.info(f"Plot saved to {save_path}")
    else:
        plt.show()
    plt.close()

def plot_beeswarm(
    df: pd.DataFrame,
    x: str,
    y: str,
    hue: Optional[str] = None,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None
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
    """
    logger = logging.getLogger(__name__)
    plt.figure(figsize=(15, 10))
    sns.swarmplot(data=df, x=x, y=y, hue=hue, legend=False)
    plt.title(title or 'Beeswarm Plot', fontsize=22)
    plt.xlabel(x, fontsize=18)
    plt.ylabel(y, fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if export and output:
        plt.savefig(output, dpi=300, format='png')
        logger.info(f"Beeswarm plot saved to {output}")
    plt.show()
    plt.close()

def plot_histogram(
    df: pd.DataFrame,
    score_col: str,
    group_col: str = 'annotate_dna',
    export: bool = False,
    output: Optional[str] = None
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
    """
    logger = logging.getLogger(__name__)
    min_score = df[score_col].min()
    max_score = df[score_col].max()
    bins = np.linspace(min_score, max_score, 51)
    unique_values = df[group_col].unique()
    colors = plt.cm.jet(np.linspace(0, 1, len(unique_values)))
    annotated_values = ['spike-in', 'wt_dna']
    other_values = [v for v in unique_values if v not in annotated_values]
    for value in other_values:
        color = colors[list(unique_values).index(value)]
        subset = df[df[group_col] == value]
        plt.hist(subset[score_col], bins=bins, color=color, alpha=0.5, label=str(value))
    for value in annotated_values:
        if value in unique_values:
            color = colors[list(unique_values).index(value)]
            subset = df[df[group_col] == value]
            plt.hist(subset[score_col], bins=bins, color=color, alpha=0.75, label=str(value), edgecolor='black')
    plt.title('Histogram of DNA level Scores')
    plt.xlabel('Activity Score')
    plt.ylabel('Frequency')
    plt.legend(title=group_col)
    if export and output:
        plt.savefig(output, dpi=1200, format='png', transparent=True)
        logger.info(f"Histogram plot saved to {output}")
    else:
        plt.show()
    plt.close()

# plot_heatmap is now imported from sortscore.visualization.heatmaps
