"""
I/O utilities for oPool DMS activity score analysis.

This module provides functions for exporting data and statistics.

Examples
--------
>>> from sortscore.analysis.io import export_dataframe, export_stats_json
"""
import logging
import pandas as pd
import json
from typing import Any, Dict

def export_dataframe(df: pd.DataFrame, filename: str) -> None:
    """
    Export a DataFrame to a CSV file.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to export.
    filename : str
        Output CSV file path.

    Returns
    -------
    None

    Examples
    --------
    >>> export_dataframe(df, 'output.csv')
    """
    logging.info(f"Exporting DataFrame to {filename}")
    df.to_csv(filename, index=False)

def export_stats_json(stats: Dict[str, Any], filename: str) -> None:
    """
    Export statistics dictionary to a JSON file.

    Parameters
    ----------
    stats : dict
        Dictionary of statistics to export.
    filename : str
        Output JSON file path.

    Returns
    -------
    None

    Examples
    --------
    >>> export_stats_json({'mean': 1.0}, 'stats.json')
    """
    logging.info(f"Exporting stats to {filename}")
    with open(filename, 'w') as f:
        json.dump(stats, f, indent=2)
