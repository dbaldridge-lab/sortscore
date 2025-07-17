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
import pandas as pd
from dataclasses import dataclass
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)

@dataclass
class ExperimentConfig:
    """
    Dataclass for experiment configuration and data loading.
    """
    experiment_name: str
    experiment_setup_file: str
    wt_seq: str
    variant_type: str
    max_pos: int
    min_pos: int
    output_dir: Optional[str] = None
    other_params: Optional[Dict[str, Any]] = None
    counts: Optional[Dict[int, Dict[int, pd.DataFrame]]] = None
    median_gfp: Optional[Dict[int, Dict[int, float]]] = None
    total_reads: Optional[Dict[int, Dict[int, int]]] = None
    
    # Analysis parameters with defaults
    bins_required: int = 1
    reps_required: int = 1  
    avg_method: str = 'rep-weighted'
    minread_threshold: int = 0
    barcoded: bool = False
    
    @property
    def num_aa(self) -> int:
        """Calculate number of amino acids from min_pos and max_pos."""
        return self.max_pos - self.min_pos + 1

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
        # Build kwargs only for fields that exist in JSON
        args = {}
        
        # Required fields
        for field in ['experiment_name', 'experiment_setup_file', 'wt_seq', 'variant_type', 'max_pos', 'min_pos']:
            if field in data:
                args[field] = data[field]
        
        # Optional fields (only add if present in JSON to preserve dataclass defaults)
        for field in ['output_dir', 'bins_required', 'reps_required', 'avg_method', 'minread_threshold', 'barcoded']:
            if field in data:
                args[field] = data[field]
        
        # Other parameters
        handled_keys = {'experiment_name', 'experiment_setup_file', 'wt_seq', 'variant_type', 'max_pos', 'min_pos',
                       'output_dir', 'bins_required', 'reps_required', 'avg_method', 'minread_threshold', 'barcoded'}
        other_params = {k: v for k, v in data.items() if k not in handled_keys}
        if other_params:
            args['other_params'] = other_params
        
        return ExperimentConfig(**args)

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
        total_reads = {}
        
        for _, row in setup_df.iterrows():
            rep = int(row['Replicate'])
            bin_ = int(row['Bin'])
            count_file = str(row['Path']).strip()
            if 'Read Counts (CSV)' in row:
                count_file = str(row['Read Counts (CSV)']).strip()
                
            gfp = float(row['Median GFP'])
            
            # Load total reads if available (for sample read depth normalization)
            if 'Read Count' in row:
                total_read_count = int(row['Read Count'])
                total_reads.setdefault(rep, {})[bin_] = total_read_count
            
            try:
                df = pd.read_csv(count_file, sep=None, engine='python', header=None)
                # Standardize column names: first column is variant sequence, second is count
                if len(df.columns) >= 2:
                    df.columns = ['variant_seq', 'count'] + list(df.columns[2:])
                else:
                    logging.error(f"Count file {count_file} does not have at least 2 columns")
                    continue
            except Exception as e:
                logging.error(f"Failed to load count file {count_file}: {e}")
                continue
            counts.setdefault(rep, {})[bin_] = df
            median_gfp.setdefault(rep, {})[bin_] = gfp
            
        self.counts = counts
        self.median_gfp = median_gfp
        
        # Set total_reads if we loaded any
        if total_reads:
            self.total_reads = total_reads
            logging.info(f"Loaded total read counts for proper normalization: {len(total_reads)} replicates")
    
    def set_total_reads(self, total_reads: Dict[int, Dict[int, int]]) -> None:
        """
        Set total sequencing reads for proper normalization (controlling for sequencing depth).
        
        Parameters
        ----------
        total_reads : dict
            Nested dict of total sequencing reads: total_reads[rep][bin] = int.
            These should be the total reads from sequencing before any filtering.
        """
        self.total_reads = total_reads

    def annotate_counts(self, wt_ref_seq: str) -> None:
        """
        Annotate all counts DataFrames with DNA sequence differences and difference counts if variant_type is 'dna'.

        Parameters
        ----------
        wt_ref_seq : str
            Wild-type DNA sequence to compare against.
        """
        from sortscore.sequence_parsing import compare_to_reference, compare_codon_lists
        if self.variant_type != 'dna' or not self.counts:
            return
        for rep_dict in self.counts.values():
            for df in rep_dict.values():
                df['dna_seq_diff'] = df['variant_seq'].apply(lambda x: compare_to_reference(wt_ref_seq, x))
                df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
                df['dna_diff_count'] = df['dna_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
                
                # Add codon diff for DNA variant types
                df['codon_diff'] = df['variant_seq'].apply(lambda x: compare_codon_lists(wt_ref_seq, x))
                df['codon_diff'] = df['codon_diff'].fillna('')
                df['codon_diff_count'] = df['codon_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)

    def annotate_variants(self, df, wt_ref_seq: str, wt_aa_seq: str, variant_type: str) -> None:
        """
        Annotate a variant count DataFrame with sequence difference columns based on variant type.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to annotate (will be modified in place).
        wt_ref_seq : str
            Wild-type DNA sequence.
        wt_aa_seq : str
            Wild-type amino acid sequence.
        variant_type : str
            'aa' for amino acid, 'dna' for nucleotide/codon. If 'dna', all annotations are added.
        """
        from sortscore.sequence_parsing import compare_to_reference, translate_dna
        
        # AA annotation (always added)
        if variant_type == 'aa':
            df['aa_seq'] = df['variant_seq']  # AA sequences provided directly
        else:  # variant_type == 'dna'
            df['aa_seq'] = df['variant_seq'].apply(translate_dna)  # Translate DNA to AA
            
        df['aa_seq_diff'] = df['aa_seq'].apply(lambda x: compare_to_reference(wt_aa_seq, x))
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        df['aa_diff_count'] = df['aa_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
        
        # DNA annotation (only for DNA variant type)
        if variant_type == 'dna':
            df['dna_seq_diff'] = df['variant_seq'].apply(lambda x: compare_to_reference(wt_ref_seq, x))
            df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
            df['dna_diff_count'] = df['dna_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
