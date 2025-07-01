"""
Parameter definitions for oPool DMS activity score analysis.

This module contains default parameters and utility functions for analysis workflows.

Examples
--------
>>> from sortscore.analysis.parameters import get_default_parameters
>>> params = get_default_parameters()
"""
from typing import Dict, Any

def get_default_parameters() -> Dict[str, Any]:
    """
    Get default analysis parameters.

    Returns
    -------
    params : dict
        Dictionary of default parameters for the analysis.

    Examples
    --------
    >>> params = get_default_parameters()
    >>> params['b']
    3
    """
    return {
        'b': 3,  # Number of bins required
        'r': 3,  # Number of replicates required
        'avg_method': 'rep-weighted',  # Averaging method
        'minread_threshold': 0,  # Minimum reads per million
    }
