"""
Utility functions for activity score analysis.

This module provides general-purpose utilities such as suffix generation and output directory creation.

Examples
--------
>>> from sortscore.analysis.utils import make_export_suffix
"""
import os
import sys
import logging
import shutil
from typing import Any, Optional, List

def make_export_suffix(experiment_name: str, b: int, minread_threshold: int, date_str: str, max_cv: Optional[float] = None) -> str:
    """
    Generate a suffix for export filenames.

    Parameters
    ----------
    experiment_name : str
        Experiment name identifier.
    b : int
        Number of bins.
    minread_threshold : int
        Minimum reads per million.
    date_str : str
        Date string (YYYYMMDD).
    max_cv : float, optional
        Maximum coefficient of variation filter.

    Returns
    -------
    suffix : str
        Export filename suffix.

    Examples
    --------
    >>> make_export_suffix('test', 3, 0, '20250701')
    'test_3-bins_0-minreads_20250701'
    >>> make_export_suffix('test', 3, 0, '20250701', 0.5)
    'test_3-bins_0-minreads_0.5-cv_20250701'
    """
    base_suffix = f'{experiment_name}_{b}-bins_{minread_threshold}-minreads'
    if max_cv is not None:
        base_suffix += f'_{max_cv}-cv'
    return f'{base_suffix}_{date_str}'

def ensure_output_subdirs(output_dir: str, subdirs: list[str] = ["scores", "figures"]) -> None:
    """
    Ensure that specified subdirectories exist within the output directory.
    
    This function handles all error cases internally and will exit the program 
    with an error message if directory creation fails.

    Parameters
    ----------
    output_dir : str
        Path to the main output directory.
    subdirs : list of str, optional
        List of subdirectory names to create (default: ["scores", "figures"]).

    Examples
    --------
    >>> ensure_output_subdirs('/path/to/output')
    # Creates /path/to/output/scores and /path/to/output/figures if they do not exist
    
    Raises
    ------
    SystemExit
        If directory creation fails due to permissions or other errors.
    """
    try:
        for subdir in subdirs:
            path = os.path.join(output_dir, subdir)
            os.makedirs(path, exist_ok=True)
    except PermissionError as e:
        logging.error(f"Permission denied: cannot create output directories in '{output_dir}'. {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Failed to create output directories in '{output_dir}': {e}")
        sys.exit(1)


def cleanup_individual_experiment_files(experiment_configs: List[str]) -> None:
    """
    Clean up individual experiment output files after batch processing.
    
    This function removes the scores and figures directories from individual
    experiment output directories, typically used after batch processing has
    combined all results into a unified output.
    
    Parameters
    ----------
    experiment_configs : List[str]
        List of paths to individual experiment config files
        
    Examples
    --------
    >>> config_files = ['exp1/config.json', 'exp2/config.json'] 
    >>> cleanup_individual_experiment_files(config_files)
    """
    from sortscore.utils.load_experiment import ExperimentConfig
    
    for config_path in experiment_configs:
        try:
            # Load experiment config to get output directory
            experiment = ExperimentConfig.from_json(config_path)
            if experiment.output_dir and os.path.exists(experiment.output_dir):
                # Remove scores and figures directories
                scores_dir = os.path.join(experiment.output_dir, 'scores')
                figures_dir = os.path.join(experiment.output_dir, 'figures')
                
                if os.path.exists(scores_dir):
                    shutil.rmtree(scores_dir)
                    logging.debug(f"Removed {scores_dir}")
                
                if os.path.exists(figures_dir):
                    shutil.rmtree(figures_dir)
                    logging.debug(f"Removed {figures_dir}")
        except Exception as e:
            logging.warning(f"Failed to cleanup files for {config_path}: {e}")
