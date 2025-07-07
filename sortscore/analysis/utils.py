"""
Utility functions for activity score analysis.

This module provides general-purpose utilities such as date formatting and suffix generation.

Examples
--------
>>> from sortscore.analysis.utils import get_current_date, make_export_suffix
>>> get_current_date()
'20250701'
"""
import datetime
import os
from typing import Any

def get_current_date() -> str:
    """
    Get the current date as a string in YYYYMMDD format.

    Returns
    -------
    date_str : str
        Current date in YYYYMMDD format.

    Examples
    --------
    >>> get_current_date()
    '20250701'
    """
    return datetime.datetime.now().strftime("%Y%m%d")

def make_export_suffix(submission: str, b: int, minread_threshold: int, date_str: str) -> str:
    """
    Generate a suffix for export filenames.

    Parameters
    ----------
    submission : str
        Submission identifier.
    b : int
        Number of bins.
    minread_threshold : int
        Minimum reads per million.
    date_str : str
        Date string (YYYYMMDD).

    Returns
    -------
    suffix : str
        Export filename suffix.

    Examples
    --------
    >>> make_export_suffix('test', 3, 0, '20250701')
    'test_3-bins_0-minreads_20250701'
    """
    return f'{submission}_{b}-bins_{minread_threshold}-minreads_{date_str}'

def ensure_output_subdirs(output_dir: str, subdirs: list[str] = ["scores", "figures"]) -> None:
    """
    Ensure that specified subdirectories exist within the output directory.

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
    """
    for subdir in subdirs:
        path = os.path.join(output_dir, subdir)
        os.makedirs(path, exist_ok=True)
