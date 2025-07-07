"""
Configuration and constants for Sort-seq variant analysis.

This module contains configuration variables and data structures for different experiment submissions.

Examples
--------
>>> from sortscore.analysis.config import get_submission_config
>>> config = get_submission_config('example_submission')
"""
import json
import logging
from dataclasses import dataclass
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

# Example structure for submission-specific configuration
def get_submission_config(submission: str) -> Optional[Dict[str, Any]]:
    """
    Retrieve configuration for a given experiment submission.

    Parameters
    ----------
    submission : str
        The submission identifier (e.g., 'example_submission').

    Returns
    -------
    config : dict or None
        Dictionary containing configuration for the submission, or None if not found.

    Examples
    --------
    >>> config = get_submission_config('example_submission')
    >>> config['gfp_rep1_bin1']
    125
    """
    configs = {
        'example_submission': {
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
            'project_dir': '/path/to/project/results/example_submission',
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

@dataclass
class ExperimentConfig:
    submission: str
    experiment_setup_file: str
    wt_seq: str
    mutant_type: str
    num_aa: int
    min_pos: int
    median_gfp: Optional[float] = None
    output_dir: Optional[str] = None
    other_params: Optional[Dict[str, Any]] = None

    @staticmethod
    def from_json(json_path: str) -> 'ExperimentConfig':
        """
        Load experiment configuration from a JSON file.

        Parameters
        ----------
        json_path : str
            Path to the JSON config file.

        Returns
        -------
        config : ExperimentConfig
            Loaded experiment configuration.
        """
        with open(json_path, 'r') as f:
            data = json.load(f)
        # Extract known fields, pass the rest to other_params
        known_fields = {k: data[k] for k in ['submission', 'experiment_setup_file', 'wt_seq', 'mutant_type', 'num_aa', 'min_pos'] if k in data}
        median_gfp = data.get('median_gfp')
        output_dir = data.get('output_dir')
        other_params = {k: v for k, v in data.items() if k not in known_fields and k not in ['median_gfp', 'output_dir']}
        return ExperimentConfig(
            **known_fields,
            median_gfp=median_gfp,
            output_dir=output_dir,
            other_params=other_params if other_params else None
        )
