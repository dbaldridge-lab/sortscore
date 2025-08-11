"""
Visualization functions and utilities for Sort-seq variant analysis.

This module provides the main plotting interface and utility functions for Sort-seq analysis.
"""

import logging
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional, List, Any, Dict, Tuple, Union
from pathlib import Path

# Import functions from specialized modules for backward compatibility
from sortscore.visualization.heatmaps import plot_heatmap
from sortscore.visualization.distributions import plot_activity_score_distribution
from sortscore.visualization.beeswarms import plot_beeswarm
from sortscore.visualization.histograms import plot_histogram

# Standard font sizes for consistent plotting
FONT_SIZES = {
    'title': 18,
    'subtitle': 16, 
    'axis_label': 14,
    'tick_label': 12,
    'legend': 12,
    'annotation': 10,
    'small': 8
}

# Standard figure sizes
FIG_SIZES = {
    'small': (8, 6),
    'medium': (12, 8),
    'large': (16, 10),
    'wide': (16, 6),
    'tall': (8, 12),
    'square': (10, 10)
}


def get_color_palette(
    palette_type: str = 'qualitative',
    n_colors: int = 10,
    custom_colors: Optional[List[str]] = None
) -> List[str]:
    """
    Get a color palette for plotting.
    
    Parameters
    ----------
    palette_type : str, default 'qualitative'
        Type of palette ('qualitative', 'sequential', 'diverging', 'pastel').
    n_colors : int, default 10
        Number of colors needed.
    custom_colors : list of str, optional
        Custom color list to use instead.
        
    Returns
    -------
    list of str
        List of color codes.
    """
    if custom_colors:
        return custom_colors[:n_colors] * (n_colors // len(custom_colors) + 1)
    
    palettes = {
        'qualitative': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
        'sequential': plt.cm.viridis(np.linspace(0, 1, n_colors)),
        'diverging': plt.cm.RdYlBu_r(np.linspace(0, 1, n_colors)),
        'pastel': ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6',
                  '#ffffcc', '#e5d8bd', '#fddaec', '#f2f2f2', '#b3b3b3']
    }
    
    if palette_type in palettes:
        colors = palettes[palette_type]
        if hasattr(colors, 'shape'):  # Matplotlib colormap array
            return ['#%02x%02x%02x' % tuple((colors[i][:3] * 255).astype(int)) 
                   for i in range(len(colors))]
        else:
            return colors[:n_colors] * (n_colors // len(colors) + 1)
    else:
        # Fallback to seaborn palette
        return sns.color_palette("husl", n_colors).as_hex()


def save_plot(
    filename: str,
    output_dir: str = '.',
    formats: List[str] = ['png'],
    dpi: int = 300,
    bbox_inches: str = 'tight',
    transparent: bool = False
) -> List[str]:
    """
    Save plot in multiple formats with consistent settings.
    
    Parameters
    ----------
    filename : str
        Base filename without extension.
    output_dir : str, default '.'
        Output directory path.
    formats : list of str, default ['png']
        List of file formats ('png', 'svg', 'pdf', 'eps').
    dpi : int, default 300
        Resolution for raster formats.
    bbox_inches : str, default 'tight'
        Bounding box setting for saving.
    transparent : bool, default False
        Whether to use transparent background.
        
    Returns
    -------
    list of str
        List of saved file paths.
    """
    logger = logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_files = []
    for fmt in formats:
        filepath = output_dir / f"{filename}.{fmt}"
        
        save_kwargs = {
            'dpi': dpi if fmt in ['png', 'jpg', 'tiff'] else None,
            'bbox_inches': bbox_inches,
            'transparent': transparent,
            'facecolor': 'white' if not transparent else 'none',
            'edgecolor': 'none'
        }
        
        # Remove None values
        save_kwargs = {k: v for k, v in save_kwargs.items() if v is not None}
        
        plt.savefig(filepath, format=fmt, **save_kwargs)
        saved_files.append(str(filepath))
        logger.info(f"Plot saved: {filepath}")
    
    return saved_files


def setup_plot_style(
    style: str = 'default',
    font_scale: float = 1.0,
    grid: bool = True,
    despine: bool = True
) -> None:
    """
    Set up consistent plot styling.
    
    Parameters
    ----------
    style : str, default 'default'
        Plot style ('default', 'publication', 'presentation').
    font_scale : float, default 1.0
        Scaling factor for font sizes.
    grid : bool, default True
        Whether to show grid lines.
    despine : bool, default True
        Whether to remove top and right spines.
    """
    if style == 'publication':
        plt.style.use('seaborn-v0_8-whitegrid' if grid else 'seaborn-v0_8-white')
        plt.rcParams.update({
            'font.size': 12 * font_scale,
            'axes.titlesize': 14 * font_scale,
            'axes.labelsize': 12 * font_scale,
            'xtick.labelsize': 10 * font_scale,
            'ytick.labelsize': 10 * font_scale,
            'legend.fontsize': 10 * font_scale,
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight'
        })
    elif style == 'presentation':
        plt.rcParams.update({
            'font.size': 16 * font_scale,
            'axes.titlesize': 20 * font_scale,
            'axes.labelsize': 18 * font_scale,
            'xtick.labelsize': 14 * font_scale,
            'ytick.labelsize': 14 * font_scale,
            'legend.fontsize': 14 * font_scale,
            'lines.linewidth': 2,
            'axes.linewidth': 1.5
        })
    
    if despine:
        sns.despine()


def format_axes(
    ax: plt.Axes,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    grid: bool = True,
    grid_alpha: float = 0.3,
    spines_to_remove: List[str] = ['top', 'right']
) -> None:
    """
    Apply consistent axis formatting.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to format.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    title : str, optional
        Plot title.
    grid : bool, default True
        Whether to show grid.
    grid_alpha : float, default 0.3
        Grid transparency.
    spines_to_remove : list of str, default ['top', 'right']
        Which spines to remove.
    """
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=FONT_SIZES['axis_label'])
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FONT_SIZES['axis_label'])
    if title:
        ax.set_title(title, fontsize=FONT_SIZES['title'])
    
    if grid:
        ax.grid(True, alpha=grid_alpha)
    
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    
    ax.tick_params(labelsize=FONT_SIZES['tick_label'])


def create_color_legend(
    colors: List[str],
    labels: List[str],
    title: Optional[str] = None,
    loc: str = 'upper right',
    bbox_to_anchor: Optional[Tuple[float, float]] = None
) -> plt.Legend:
    """
    Create a custom color legend.
    
    Parameters
    ----------
    colors : list of str
        List of colors.
    labels : list of str
        List of labels.
    title : str, optional
        Legend title.
    loc : str, default 'upper right'
        Legend location.
    bbox_to_anchor : tuple, optional
        Bbox anchor coordinates.
        
    Returns
    -------
    matplotlib.legend.Legend
        Legend object.
    """
    legend_elements = [patches.Patch(facecolor=color, label=label)
                      for color, label in zip(colors, labels)]
    
    legend_kwargs = {
        'handles': legend_elements,
        'loc': loc,
        'fontsize': FONT_SIZES['legend']
    }
    
    if title:
        legend_kwargs['title'] = title
    
    if bbox_to_anchor:
        legend_kwargs['bbox_to_anchor'] = bbox_to_anchor
    
    return plt.legend(**legend_kwargs)


def create_subplot_grid(
    nrows: int,
    ncols: int,
    figsize: Optional[Tuple[float, float]] = None,
    subplot_titles: Optional[List[str]] = None,
    main_title: Optional[str] = None,
    share_x: bool = False,
    share_y: bool = False
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create a grid of subplots with consistent formatting.
    
    Parameters
    ----------
    nrows : int
        Number of subplot rows.
    ncols : int
        Number of subplot columns.
    figsize : tuple, optional
        Figure size (width, height).
    subplot_titles : list of str, optional
        Titles for each subplot.
    main_title : str, optional
        Main figure title.
    share_x : bool, default False
        Whether to share x-axes.
    share_y : bool, default False
        Whether to share y-axes.
        
    Returns
    -------
    tuple
        Figure and axes array.
    """
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                           sharex=share_x, sharey=share_y)
    
    # Ensure axes is always an array
    if nrows * ncols == 1:
        axes = np.array([axes])
    elif nrows == 1 or ncols == 1:
        axes = axes.reshape(nrows, ncols)
    
    # Set subplot titles
    if subplot_titles:
        for i, (ax, title) in enumerate(zip(axes.flat, subplot_titles)):
            ax.set_title(title, fontsize=FONT_SIZES['subtitle'])
    
    # Set main title
    if main_title:
        fig.suptitle(main_title, fontsize=FONT_SIZES['title'])
    
    plt.tight_layout()
    return fig, axes