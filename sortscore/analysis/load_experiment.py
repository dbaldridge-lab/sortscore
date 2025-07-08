"""
Configuration and constants for Sort-seq variant analysis.

This module provides the ExperimentConfig dataclass for experiment configuration and data loading.

Examples
--------
>>> from sortscore.analysis.config import ExperimentConfig
>>> config = ExperimentConfig.from_json('config.json')
>>> config.load_counts()
>>> df = config.counts[1][2]  # Access replicate 1, bin 2
"""
import json
import logging
from dataclasses import dataclass
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

@dataclass
class ExperimentConfig:
    """
    Dataclass for experiment configuration and data loading.
    """
    submission: str
    experiment_setup_file: str
    wt_seq: str
    mutant_type: str
    num_aa: int
    min_pos: int
    median_gfp: Optional[float] = None
    output_dir: Optional[str] = None
    other_params: Optional[Dict[str, Any]] = None
    counts: Optional[Dict[int, Dict[int, 'pd.DataFrame']]] = None

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

    def load_counts(self) -> None:
        """
        Load count DataFrames for all replicates and bins as specified in the experiment setup file.

        This populates the `counts` attribute as a nested dictionary: counts[rep][bin] = DataFrame.

        Returns
        -------
        None
        """
        import pandas as pd
        try:
            setup_df = pd.read_csv(self.experiment_setup_file)
        except Exception as e:
            raise RuntimeError(f"Failed to load experiment setup file: {e}")
        counts = {}
        for _, row in setup_df.iterrows():
            rep = int(row['replicate'])
            bin_ = int(row['bin'])
            count_file = row['count_file']
            try:
                df = pd.read_csv(count_file, sep=None, engine='python')
            except Exception as e:
                logging.error(f"Failed to load count file {count_file}: {e}")
                continue
            counts.setdefault(rep, {})[bin_] = df
        self.counts = counts
