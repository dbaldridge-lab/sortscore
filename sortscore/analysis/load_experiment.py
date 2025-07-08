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
    output_dir: Optional[str] = None
    other_params: Optional[Dict[str, Any]] = None
    counts: Optional[Dict[int, Dict[int, 'pd.DataFrame']]] = None
    median_gfp: Optional[Dict[int, Dict[int, float]]] = None

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
        output_dir = data.get('output_dir')
        other_params = {k: v for k, v in data.items() if k not in known_fields and k != 'output_dir'}
        return ExperimentConfig(
            **known_fields,
            output_dir=output_dir,
            other_params=other_params if other_params else None
        )

    def load_counts(self) -> None:
        """
        Load count DataFrames for all replicates and bins as specified in the experiment setup file.

        This populates the `counts` attribute as a nested dictionary: counts[rep][bin] = DataFrame.
        Also populates the `median_gfp` attribute as a nested dictionary: median_gfp[rep][bin] = float.

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
        median_gfp = {}
        for _, row in setup_df.iterrows():
            rep = int(row['Replicate'])
            bin_ = int(row['Bin'])
            count_file = row['Read Counts (CSV)']
            gfp = float(row['Median GFP'])
            try:
                df = pd.read_csv(count_file, sep=None, engine='python')
            except Exception as e:
                logging.error(f"Failed to load count file {count_file}: {e}")
                continue
            counts.setdefault(rep, {})[bin_] = df
            median_gfp.setdefault(rep, {})[bin_] = gfp
        self.counts = counts
        self.median_gfp = median_gfp

    def annotate_counts(self, wt_ref_seq: str) -> None:
        """
        Annotate all counts DataFrames with DNA sequence differences and difference counts if mutant_type is 'dna'.

        Parameters
        ----------
        wt_ref_seq : str
            Wild-type DNA sequence to compare against.
        """
        from sortscore.sequence_parsing import compare_to_reference
        if self.mutant_type != 'dna' or not self.counts:
            return
        for rep_dict in self.counts.values():
            for df in rep_dict.values():
                df['dna_seq_diff'] = df['seq'].apply(lambda x: compare_to_reference(wt_ref_seq, x))
                df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
                df['dna_diff_count'] = df['dna_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)

    def annotate_variants(self, df, wt_ref_seq: str, wt_aa_seq: str, mutant_type: str) -> None:
        """
        Annotate a variant count DataFrame with sequence difference columns based on mutant type.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to annotate (will be modified in place).
        wt_ref_seq : str
            Wild-type DNA sequence.
        wt_aa_seq : str
            Wild-type amino acid sequence.
        mutant_type : str
            'aa' for amino acid, 'dna' for nucleotide/codon. If 'dna', all annotations are added.
        """
        from sortscore.sequence_parsing import compare_to_reference, translate_dna
        # AA annotation (always added for aa or dna)
        df['aa_seq'] = df['variant_seq'].apply(translate_dna)
        df['aa_seq_diff'] = df['aa_seq'].apply(lambda x: compare_to_reference(wt_aa_seq, x))
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        df['aa_diff_count'] = df['aa_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
        if mutant_type == 'dna':
            # DNA annotation
            df['dna_seq_diff'] = df['variant_seq'].apply(lambda x: compare_to_reference(wt_ref_seq, x))
            df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
            df['dna_diff_count'] = df['dna_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
