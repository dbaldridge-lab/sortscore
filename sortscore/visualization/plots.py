"""
Visualization functions for Sort-seq variant analysis.

This module provides plotting functions for variant activity score data, including beeswarm, heatmap, and histogram plots.

Examples
--------
>>> from sortscore.visualization.plots import plot_activity_score_distribution, plot_beeswarm, plot_heatmap, plot_histogram
"""
import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, List, Any
from sortscore.visualization.heatmap_matrix import make_dms_matrix, fill_wt, make_col_avg_df, get_dropout
from sortscore.analysis.load_experiment import ExperimentConfig


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

def plot_heatmap(
    data: pd.DataFrame,
    score_col: str,
    experiment: ExperimentConfig,
    wt_score: float = 1.0,
    fig_size: str = 'small',
    export: bool = False,
    output: Optional[str] = None,
    tick_values: Optional[List[float]] = None,
    tick_labels: Optional[List[str]] = None,
    motif_indices: Optional[List[int]] = None,
    row_avg: bool = False,
    title: Optional[str] = None
) -> None:
    """
    Plot a DMS heatmap using a matrix of activity scores.

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
    """
    logger = logging.getLogger(__name__)
    dms_matrix = make_dms_matrix(
        data,
        score_col,
        experiment.num_aa,
        experiment.wt_seq,
        experiment.min_pos,
        experiment.variant_type
    )
    dropout_num, dropout_percent = get_dropout(dms_matrix)
    heatmap_df = fill_wt(dms_matrix, wt_score)
    col_avg_df = make_col_avg_df(heatmap_df)
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.gridspec import GridSpec
    nan_mask = dms_matrix.isnull()
    wt_mask = dms_matrix == 'WT'
    if fig_size == 'small':
        fig = plt.figure(figsize=(16.5, 12))
        tick_freq = 2
    elif fig_size == 'large':
        fig = plt.figure(figsize=(30, 10))
        tick_freq = 5
    elif fig_size == 'long':
        fig = plt.figure(figsize=(30, 25))
        tick_freq = 5
    if row_avg:
        row_avg_df = pd.DataFrame(heatmap_df.mean(axis=1), columns=['Avg'])
        gs = GridSpec(2,3, width_ratios=[1, 35, 1], height_ratios=[1, 45], hspace=0.03, wspace=0.03)
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[1, 1])
        cax = fig.add_subplot(gs[1, 2])
        ax3 = fig.add_subplot(gs[1, 0])
    else:
        gs = GridSpec(2, 2, width_ratios=[35, 1], height_ratios=[1, 20], hspace=0.03, wspace=0.03)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        cax = fig.add_subplot(gs[1, 1])
    min_val = heatmap_df.min().min()
    max_val = heatmap_df.max().max()
    norm = plt.Normalize(vmin=min_val, vmax=max_val)
    cmap = plt.cm.magma
    sns.heatmap(col_avg_df, annot=False, cmap=cmap, cbar=False, ax=ax1, norm=norm)
    ax = sns.heatmap(heatmap_df, cmap=cmap, cbar=False, ax=ax2)
    if nan_mask.any().any():
        sns.heatmap(nan_mask, annot=False, cmap=cmap, cbar=False, mask=~nan_mask, ax=ax2)
    if row_avg:
        sns.heatmap(row_avg_df, cmap=cmap, cbar=False, ax=ax3, norm=norm)
    for i in range(len(nan_mask)):
        for j in range(len(nan_mask.columns)):
            if nan_mask.iloc[i, j]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, hatch='//', edgecolor='lightgray', facecolor='white'))
    wt_indices = np.where(wt_mask)
    ax2.scatter(wt_indices[1] + 0.5, wt_indices[0] + 0.5, color='white', s=30, alpha=0.5)
    if motif_indices:
        pass
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    if tick_values is not None and tick_labels is not None:
        cbar = plt.colorbar(sm, cax=cax, ticks=tick_values)
        cbar.ax.set_yticklabels(tick_labels, fontsize=18)
    else:
        cbar = plt.colorbar(sm, cax=cax)
    x_labels = [str(i) for i in range(experiment.min_pos, experiment.min_pos+experiment.num_aa)]
    ax.set_xticks([i + 0.5 for i in range(0, len(x_labels), tick_freq)])
    ax.set_xticklabels([x_labels[i] for i in range(0, len(x_labels), tick_freq)], rotation=0)
    plot_title = title or f'DMS Heatmap - Dropout {dropout_num} variant ({dropout_percent}%)'
    if fig_size == 'small':
        fig.suptitle(plot_title, fontsize=28, y=0.89)
    elif fig_size == 'large':
        fig.suptitle(plot_title, fontsize=30, y=0.91)
    fig.subplots_adjust(top=0.85, bottom=0.15)
    ax1.set_ylabel('Avg', fontsize=20, labelpad=8, ha='center')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_xlabel('Residue Sequence Number', fontsize=24)
    ax2.set_ylabel('Variant Amino Acid', fontsize=24)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    if row_avg:
        ax3.set_ylabel('Variant Amino Acid', fontsize=24, labelpad=8, ha='center')
        ax3.set_yticks(ax2.get_yticks())
        ax3.set_yticklabels(ax2.get_yticklabels())
        ax3.tick_params(axis='y', which='major', labelsize=20)
        ax2.set_ylabel('')
        ax2.set_yticks([])
        ax3.set_xticks([])
    if export and output:
        plt.savefig(output, dpi=300, format='png', transparent=True)
        logger.info(f"Heatmap plot saved to {output}")
    else:
        plt.show()
    plt.close()

# Note: Heatmap and matrix utilities are complex and depend on sequence parsing and stats modules.
# For a full refactor, these should be moved to a separate module and depend on the new package structure.
