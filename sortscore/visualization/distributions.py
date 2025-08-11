"""
Activity score distribution analysis and visualization functions for Sort-seq variant analysis.

This module provides functions for analyzing and visualizing activity score distributions,
including histograms, density plots, Q-Q plots, and box plots.

Examples
--------
>>> from sortscore.visualization.distributions import plot_activity_score_distribution, plot_density_comparison
>>> plot_activity_score_distribution(df, score_col='avgscore_rep_weighted')
"""

import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
from typing import Optional, List, Dict, Any
from matplotlib.patches import Rectangle

def plot_activity_score_distribution(
    df: pd.DataFrame,
    score_col: str = 'avgscore',
    bins: int = 50,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show_stats: bool = True,
    overlay_normal: bool = False
) -> Dict[str, float]:
    """
    Plot the distribution of activity scores with optional statistics and normal overlay.

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
    show_stats : bool, default True
        If True, display distribution statistics on plot.
    overlay_normal : bool, default False
        If True, overlay fitted normal distribution.

    Returns
    -------
    dict
        Dictionary containing distribution statistics.
    """
    logger = logging.getLogger(__name__)
    if score_col not in df.columns:
        logger.error(f"Score column '{score_col}' not found in DataFrame.")
        raise ValueError(f"Score column '{score_col}' not found in DataFrame.")
    
    # Clean data
    clean_scores = df[score_col].dropna()
    
    # Calculate statistics
    stats_dict = {
        'mean': clean_scores.mean(),
        'median': clean_scores.median(),
        'std': clean_scores.std(),
        'skew': stats.skew(clean_scores),
        'kurtosis': stats.kurtosis(clean_scores),
        'n': len(clean_scores)
    }
    
    plt.figure(figsize=(10, 6))
    
    # Create histogram
    n_values, bin_edges, patches = plt.hist(clean_scores, bins=bins, density=True, 
                                          alpha=0.7, color='skyblue', edgecolor='black')
    
    # Overlay normal distribution if requested
    if overlay_normal:
        x = np.linspace(clean_scores.min(), clean_scores.max(), 100)
        normal_pdf = stats.norm.pdf(x, stats_dict['mean'], stats_dict['std'])
        plt.plot(x, normal_pdf, 'r-', linewidth=2, label='Normal fit')
        plt.legend()
    
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.title(title or f'Activity Score Distribution ({score_col})', fontsize=16)
    plt.grid(True, alpha=0.3)
    
    # Add statistics text
    if show_stats:
        stats_text = (f"Mean: {stats_dict['mean']:.3f}\n"
                     f"Median: {stats_dict['median']:.3f}\n"
                     f"Std: {stats_dict['std']:.3f}\n"
                     f"Skew: {stats_dict['skew']:.3f}\n"
                     f"n = {stats_dict['n']}")
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
                verticalalignment='top', fontsize=12,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Distribution plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()
    return stats_dict


def plot_density_comparison(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    palette: Optional[str] = 'Set2'
) -> None:
    """
    Plot kernel density estimates for different groups.
    
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
    save_path : str, optional
        If provided, save the plot to this file path.
    palette : str, optional
        Color palette for groups.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 6))
    
    # Create density plots for each group
    groups = df[group_col].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))
    
    for i, group in enumerate(groups):
        group_data = df[df[group_col] == group][score_col].dropna()
        if len(group_data) > 1:
            density = stats.gaussian_kde(group_data)
            x = np.linspace(df[score_col].min(), df[score_col].max(), 200)
            plt.plot(x, density(x), color=colors[i], linewidth=2, label=f'{group} (n={len(group_data)})')
    
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.title(title or f'Activity Score Density by {group_col}', fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Density comparison saved to {save_path}")
    else:
        plt.show()
    
    plt.close()


def plot_qq_normal(
    df: pd.DataFrame,
    score_col: str,
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> Dict[str, float]:
    """
    Create Q-Q plot against normal distribution to test normality.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_col : str
        Column name for the activity score.
    title : str, optional
        Plot title.
    save_path : str, optional
        If provided, save the plot to this file path.
        
    Returns
    -------
    dict
        Normality test results.
    """
    logger = logging.getLogger(__name__)
    
    clean_scores = df[score_col].dropna()
    
    # Perform normality tests
    shapiro_stat, shapiro_p = stats.shapiro(clean_scores[:5000])  # Limit for Shapiro-Wilk
    ks_stat, ks_p = stats.kstest(clean_scores, 'norm', 
                                args=(clean_scores.mean(), clean_scores.std()))
    
    test_results = {
        'shapiro_stat': shapiro_stat,
        'shapiro_p': shapiro_p,
        'ks_stat': ks_stat,
        'ks_p': ks_p
    }
    
    plt.figure(figsize=(8, 8))
    
    # Create Q-Q plot
    stats.probplot(clean_scores, dist="norm", plot=plt)
    plt.title(title or f'Q-Q Plot: {score_col} vs Normal', fontsize=16)
    plt.grid(True, alpha=0.3)
    
    # Add test results
    test_text = (f"Shapiro-Wilk: p = {shapiro_p:.2e}\n"
                f"Kolmogorov-Smirnov: p = {ks_p:.2e}")
    plt.text(0.02, 0.98, test_text, transform=plt.gca().transAxes,
            verticalalignment='top', fontsize=12,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Q-Q plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()
    return test_results


def plot_box_comparison(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show_outliers: bool = True
) -> None:
    """
    Create box plots comparing score distributions across groups.
    
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
    save_path : str, optional
        If provided, save the plot to this file path.
    show_outliers : bool, default True
        If True, show outlier points.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 8))
    
    # Create box plot
    box_data = []
    labels = []
    
    for group in sorted(df[group_col].unique()):
        group_scores = df[df[group_col] == group][score_col].dropna()
        if len(group_scores) > 0:
            box_data.append(group_scores)
            labels.append(f"{group}\n(n={len(group_scores)})")
    
    bp = plt.boxplot(box_data, labels=labels, patch_artist=True, 
                    showfliers=show_outliers, notch=True)
    
    # Color boxes
    colors = plt.cm.Set3(np.linspace(0, 1, len(box_data)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    plt.ylabel('Activity Score', fontsize=14)
    plt.xlabel(group_col, fontsize=14)
    plt.title(title or f'Activity Score Distribution by {group_col}', fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Box plot comparison saved to {save_path}")
    else:
        plt.show()
    
    plt.close()


def plot_multi_histogram(
    df: pd.DataFrame,
    score_cols: List[str],
    bins: int = 30,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    alpha: float = 0.7
) -> None:
    """
    Plot overlapping histograms for multiple score columns.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    score_cols : list of str
        List of column names for scores to compare.
    bins : int, default 30
        Number of bins for histograms.
    title : str, optional
        Plot title.
    save_path : str, optional
        If provided, save the plot to this file path.
    alpha : float, default 0.7
        Transparency for overlapping histograms.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 8))
    
    # Determine common bin range
    all_scores = pd.concat([df[col].dropna() for col in score_cols])
    bin_range = (all_scores.min(), all_scores.max())
    
    colors = plt.cm.Set2(np.linspace(0, 1, len(score_cols)))
    
    for i, col in enumerate(score_cols):
        clean_scores = df[col].dropna()
        plt.hist(clean_scores, bins=bins, range=bin_range, alpha=alpha,
                color=colors[i], label=f"{col} (n={len(clean_scores)})",
                density=True, edgecolor='black', linewidth=0.5)
    
    plt.xlabel('Activity Score', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.title(title or 'Activity Score Distribution Comparison', fontsize=16)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Multi-histogram saved to {save_path}")
    else:
        plt.show()
    
    plt.close()


def plot_violin_comparison(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show_points: bool = False
) -> None:
    """
    Create violin plots comparing score distributions across groups.
    
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
    save_path : str, optional
        If provided, save the plot to this file path.
    show_points : bool, default False
        If True, overlay data points.
    """
    logger = logging.getLogger(__name__)
    
    plt.figure(figsize=(12, 8))
    
    # Create violin plot
    sns.violinplot(data=df, x=group_col, y=score_col, inner='quartile')
    
    if show_points:
        sns.stripplot(data=df, x=group_col, y=score_col, size=2, 
                     color='black', alpha=0.3)
    
    plt.ylabel('Activity Score', fontsize=14)
    plt.xlabel(group_col, fontsize=14)
    plt.title(title or f'Activity Score Distribution by {group_col}', fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Violin plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()