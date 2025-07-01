"""
Configuration and constants for oPool DMS activity score analysis.

This module contains configuration variables and data structures for different oPool submissions.

Examples
--------
>>> from sortscore.analysis.config import get_submission_config
>>> config = get_submission_config('030325_oPool5b_GTAC')
"""
import logging
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)

# Example structure for submission-specific configuration
def get_submission_config(submission: str) -> Optional[Dict[str, Any]]:
    """
    Retrieve configuration for a given oPool submission.

    Parameters
    ----------
    submission : str
        The submission identifier (e.g., '030325_oPool5b_GTAC').

    Returns
    -------
    config : dict or None
        Dictionary containing configuration for the submission, or None if not found.

    Examples
    --------
    >>> config = get_submission_config('030325_oPool5b_GTAC')
    >>> config['gfp_rep1_bin1']
    125
    """
    configs = {
        '030325_oPool5b_GTAC': {
            'gfp_rep1_bin1': 125,
            'gfp_rep1_bin2': 242,
            'gfp_rep1_bin3': 953,
            'gfp_rep1_bin4': 3077,
            'gfp_rep2_bin1': 111,
            'gfp_rep2_bin2': 226,
            'gfp_rep2_bin3': 827,
            'gfp_rep2_bin4': 3282,
            'gfp_rep3_bin1': 86.9,
            'gfp_rep3_bin2': 204,
            'gfp_rep3_bin3': 711,
            'gfp_rep3_bin4': 3038,
            'read_count': [93335051 + 38013168, 102606862 + 45307832, 74357968 + 32724818, 75393971 + 33106487,
                           93094462 + 39259132, 64549755 + 30382114, 89303873 + 34775693, 100813592 + 43668775,
                           79373783 + 34415035, 93444866 + 37572151, 66659708 + 32502674, 86346094 + 36977168,
                           85468879 + 33582490, 76316244 + 33633147],
            'project_dir': '/scratch/dblab/opool/data/results/oPool5/030325_oPool5b_GTAC_allowlist-match',
            'tsv_file': [
                f'S{i}_matched_allowlist_substring_unique_counts.tsv.gz' for i in range(14)
            ],
        },
        # ... Add other submissions as needed ...
    }
    config = configs.get(submission)
    if config is None:
        logger.warning(f"No configuration found for submission: {submission}")
    return config
