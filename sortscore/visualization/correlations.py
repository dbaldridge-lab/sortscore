"""
Replicate correlation analysis and visualization functions for Sort-seq variant analysis.

This module provides functions for analyzing and visualizing correlations between experimental replicates,
including scatter plots, correlation matrices, and consistency metrics.

Examples
--------
>>> from sortscore.visualization.correlations import plot_replicate_correlation, plot_correlation_matrix
>>> plot_replicate_correlation(data, 'rep1_score', 'rep2_score')
"""

import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
from typing import Optional, List, Tuple, Dict
from matplotlib.patches import Rectangle

def plot_replicate_correlation(
    df: pd.DataFrame,
    rep1_col: str,
    rep2_col: str,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    format: str = 'png',
    dpi: int = 300,
    show_stats: bool = True
) -> Tuple[float, float]:
    """
    Plot correlation between two replicates with scatter plot and statistics.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing replicate data.
    rep1_col : str
        Column name for first replicate scores.
    rep2_col : str
        Column name for second replicate scores.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    format : str, default 'png'
        Output format ('png' or 'svg').
    dpi : int, default 300
        Resolution for saved plot.
    show_stats : bool, default True
        If True, display correlation statistics on plot.
        
    Returns
    -------
    tuple of float
        Pearson correlation coefficient and p-value.
    """
    logger = logging.getLogger(__name__)
    
    # Remove NaN values for correlation calculation
    clean_data = df[[rep1_col, rep2_col]].dropna()
    
    if len(clean_data) < 3:
        logger.warning(f"Insufficient data for correlation: {len(clean_data)} valid points")
        return np.nan, np.nan
    
    # Calculate correlation
    r, p = stats.pearsonr(clean_data[rep1_col], clean_data[rep2_col])
    
    # Create plot
    plt.figure(figsize=(8, 8))
    plt.scatter(clean_data[rep1_col], clean_data[rep2_col], alpha=0.6, s=20)
    
    # Add regression line
    z = np.polyfit(clean_data[rep1_col], clean_data[rep2_col], 1)
    p_line = np.poly1d(z)
    x_line = np.linspace(clean_data[rep1_col].min(), clean_data[rep1_col].max(), 100)
    plt.plot(x_line, p_line(x_line), "r--", alpha=0.8, linewidth=2)
    
    # Add unity line
    min_val = min(clean_data[rep1_col].min(), clean_data[rep2_col].min())
    max_val = max(clean_data[rep1_col].max(), clean_data[rep2_col].max())
    plt.plot([min_val, max_val], [min_val, max_val], 'k-', alpha=0.5, linewidth=1)
    
    plt.xlabel(rep1_col, fontsize=16)
    plt.ylabel(rep2_col, fontsize=16)
    plt.title(title or f'Replicate Correlation: {rep1_col} vs {rep2_col}', fontsize=18)
    
    # Add statistics text
    if show_stats:
        stats_text = f'r = {r:.3f}\np = {p:.2e}\nn = {len(clean_data)}'
        plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
                verticalalignment='top', fontsize=14, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=dpi, format=format, facecolor='white', edgecolor='none')
        logger.info(f"Correlation plot saved to {output}")
    else:
        plt.show()
    
    plt.close()
    return r, p


def plot_correlation_matrix(
    df: pd.DataFrame,
    score_cols: List[str],
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None,
    format: str = 'png',
    dpi: int = 300,
    method: str = 'pearson'
) -> pd.DataFrame:
    """
    Plot correlation matrix heatmap for multiple replicates.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing replicate data.
    score_cols : list of str
        List of column names for replicate scores.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
    format : str, default 'png'
        Output format ('png' or 'svg').
    dpi : int, default 300
        Resolution for saved plot.
    method : str, default 'pearson'
        Correlation method ('pearson', 'spearman', 'kendall').
        
    Returns
    -------
    pandas.DataFrame
        Correlation matrix.
    """
    logger = logging.getLogger(__name__)
    
    # Calculate correlation matrix
    corr_matrix = df[score_cols].corr(method=method)
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Create heatmap
    sns.heatmap(corr_matrix, annot=True, cmap='RdYlBu_r', center=0, 
                square=True, fmt='.3f', cbar_kws={'shrink': 0.8},
                annot_kws={'size': 12})
    
    plt.title(title or f'{method.capitalize()} Correlation Matrix', fontsize=18)
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=dpi, format=format, facecolor='white', edgecolor='none')
        logger.info(f"Correlation matrix saved to {output}")
    else:
        plt.show()
    
    plt.close()
    return corr_matrix


def plot_replicate_consistency(
    df: pd.DataFrame,
    score_cols: List[str],
    threshold: float = 0.1,
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None
) -> Dict[str, float]:
    """
    Plot replicate consistency analysis showing coefficient of variation.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing replicate data.
    score_cols : list of str
        List of column names for replicate scores.
    threshold : float, default 0.1
        CV threshold for highlighting inconsistent variants.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
        
    Returns
    -------
    dict
        Summary statistics of replicate consistency.
    """
    logger = logging.getLogger(__name__)
    
    # Calculate coefficient of variation across replicates
    replicate_data = df[score_cols].dropna()
    cv_values = replicate_data.std(axis=1) / replicate_data.mean(axis=1).abs()
    cv_values = cv_values.replace([np.inf, -np.inf], np.nan).dropna()
    
    # Summary statistics
    stats = {
        'mean_cv': cv_values.mean(),
        'median_cv': cv_values.median(),
        'high_cv_fraction': (cv_values > threshold).mean(),
        'n_variants': len(cv_values)
    }
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # CV histogram
    ax1.hist(cv_values, bins=50, alpha=0.7, edgecolor='black')
    ax1.axvline(threshold, color='red', linestyle='--', linewidth=2, 
                label=f'Threshold = {threshold}')
    ax1.set_xlabel('Coefficient of Variation', fontsize=14)
    ax1.set_ylabel('Count', fontsize=14)
    ax1.set_title('Replicate CV Distribution', fontsize=16)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # CV vs mean score
    mean_scores = replicate_data.mean(axis=1)
    ax2.scatter(mean_scores, cv_values, alpha=0.6, s=20)
    ax2.axhline(threshold, color='red', linestyle='--', linewidth=2)
    ax2.set_xlabel('Mean Score', fontsize=14)
    ax2.set_ylabel('Coefficient of Variation', fontsize=14)
    ax2.set_title('CV vs Mean Score', fontsize=16)
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle(title or 'Replicate Consistency Analysis', fontsize=18)
    plt.tight_layout()
    
    # Add statistics text
    stats_text = (f"Mean CV: {stats['mean_cv']:.3f}\n"
                  f"Median CV: {stats['median_cv']:.3f}\n"
                  f"High CV fraction: {stats['high_cv_fraction']:.1%}")
    fig.text(0.02, 0.02, stats_text, fontsize=12, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Consistency analysis saved to {output}")
    else:
        plt.show()
    
    plt.close()
    return stats


def plot_pairwise_comparisons(
    df: pd.DataFrame,
    score_cols: List[str],
    title: Optional[str] = None,
    export: bool = False,
    output: Optional[str] = None
) -> pd.DataFrame:
    """
    Create pairwise correlation plots for all replicate combinations.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing replicate data.
    score_cols : list of str
        List of column names for replicate scores.
    title : str, optional
        Plot title.
    export : bool, default False
        If True, save the plot to file.
    output : str, optional
        Output file path if export is True.
        
    Returns
    -------
    pandas.DataFrame
        Pairwise correlation statistics.
    """
    logger = logging.getLogger(__name__)
    
    n_reps = len(score_cols)
    n_pairs = n_reps * (n_reps - 1) // 2
    
    if n_pairs == 0:
        logger.warning("Need at least 2 replicates for pairwise comparison")
        return pd.DataFrame()
    
    # Calculate grid dimensions
    cols = int(np.ceil(np.sqrt(n_pairs)))
    rows = int(np.ceil(n_pairs / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 5*rows))
    if n_pairs == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    correlations = []
    pair_idx = 0
    
    for i in range(n_reps):
        for j in range(i+1, n_reps):
            rep1, rep2 = score_cols[i], score_cols[j]
            clean_data = df[[rep1, rep2]].dropna()
            
            if len(clean_data) < 3:
                continue
            
            r, p = stats.pearsonr(clean_data[rep1], clean_data[rep2])
            correlations.append({
                'rep1': rep1, 'rep2': rep2, 'correlation': r, 'p_value': p, 'n': len(clean_data)
            })
            
            ax = axes[pair_idx]
            ax.scatter(clean_data[rep1], clean_data[rep2], alpha=0.6, s=15)
            
            # Add regression line
            z = np.polyfit(clean_data[rep1], clean_data[rep2], 1)
            p_line = np.poly1d(z)
            x_line = np.linspace(clean_data[rep1].min(), clean_data[rep1].max(), 100)
            ax.plot(x_line, p_line(x_line), "r--", alpha=0.8, linewidth=2)
            
            ax.set_xlabel(rep1, fontsize=12)
            ax.set_ylabel(rep2, fontsize=12)
            ax.set_title(f'r = {r:.3f}', fontsize=14)
            ax.grid(True, alpha=0.3)
            
            pair_idx += 1
    
    # Hide unused subplots
    for idx in range(pair_idx, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(title or 'Pairwise Replicate Correlations', fontsize=18)
    plt.tight_layout()
    
    if export and output:
        plt.savefig(output, dpi=300, format='png', facecolor='white', edgecolor='none')
        logger.info(f"Pairwise correlations saved to {output}")
    else:
        plt.show()
    
    plt.close()
    
    return pd.DataFrame(correlations)