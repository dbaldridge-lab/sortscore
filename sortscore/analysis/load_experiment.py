"""
Configuration and constants for Sort-seq variant analysis.

This module provides the ExperimentConfig dataclass for experiment configuration and data loading.
It includes automatic detection and conversion of pre-annotated amino acid changes.

Supported Variant Formats
-------------------------
1. **Full Sequences**: DNA or protein sequences to be compared against the wild-type sequence.
   - DNA: Full nucleotide sequences for variant_type="dna"
   - AA: Full amino acid sequences for variant_type="aa"

TODO: Add HGVS c. notation support
2. **Pre-annotated Variants** (Auto-detected):
   - Single-letter codes: "M1M", "R98C", "P171X"
   - Three-letter codes: "Met1Met", "Arg98Cys", "Pro171Ter"
   - HGVS p. notation: "p.M1M", "p.Arg98Cys", "p.Pro171Ter"
   - Other delimiters: "M.1.M", "R-98-C", "P_171_X"

Variants are converted to a consistent internal annotation format for downstream analysis.

Examples
--------
>>> from sortscore.analysis.load_experiment import ExperimentConfig
>>> config = ExperimentConfig.from_json('config.json')
>>> config.load_counts()  # Automatically detects and converts pre-annotated formats
>>> df = config.counts[1][2]  # Access replicate 1, bin 2 read counts
"""
import json
import logging
import pandas as pd
import re
from dataclasses import dataclass
from typing import Dict, Any, Optional
from sortscore.sequence_parsing import get_reference_sequence

logger = logging.getLogger(__name__)

def is_aa_change_format(variant_seq: str) -> bool:
    """
    Detect if a variant sequence is in pre-annotated AA change format.
    
    Supports multiple formats:
    - Single-letter: "M1M", "R98C", "P171X"
    - Three-letter: "Met1Met", "Arg98Cys", "Pro171Ter"
    - HGVS p. notation: "p.M1M", "p.Arg98Cys", "p.Pro171Ter"
    
    Parameters
    ----------
    variant_seq : str
        Variant sequence to check
        
    Returns
    -------
    bool
        True if the sequence matches AA change format
    """
    # Clean the sequence - remove p. prefix, whitespace, and common separators
    clean_seq = re.sub(r'^p\.', '', variant_seq.strip())
    clean_seq = re.sub(r'[\s\.\-_]', '', clean_seq)
    
    # Three-letter AA codes
    aa_codes = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 
               'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 
               'Ter', 'Amb', 'Sec', 'Pyl']
    
    # Try three-letter format first
    aa_pattern = '|'.join(aa_codes)
    pattern_3letter = rf'^({aa_pattern})([0-9]+)({aa_pattern})$'
    if re.match(pattern_3letter, clean_seq, re.IGNORECASE):
        return True
    
    # Try single-letter format
    pattern_1letter = r'^[A-Z*X][0-9]+[A-Z*X]$'
    if re.match(pattern_1letter, clean_seq, re.IGNORECASE):
        return True
    
    return False

def parse_aa_change(variant_seq: str) -> tuple:
    """
    Parse a pre-annotated AA change into components.
    
    Supports multiple formats:
    - Single-letter: "M1M", "R98C", "P171X"
    - Three-letter: "Met1Met", "Arg98Cys", "Pro171Ter"
    - HGVS p. notation: "p.M1M", "p.Arg98Cys", "p.Pro171Ter"
    
    Parameters
    ----------
    variant_seq : str
        Variant sequence in various formats
        
    Returns
    -------
    tuple
        (ref_aa, position, alt_aa) where position is 1-based and AAs are single-letter codes
    """
    # Clean the sequence - remove p. prefix, whitespace, and common separators
    clean_seq = re.sub(r'^p\.', '', variant_seq.strip())
    clean_seq = re.sub(r'[\s\.\-_]', '', clean_seq)
    
    # Mapping from three-letter to single-letter codes
    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
        'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
        'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
        'Tyr': 'Y', 'Val': 'V', 'Ter': '*', 'Amb': 'X', 'Sec': 'U', 'Pyl': 'O'
    }
    
    # Three-letter AA codes
    aa_codes = list(aa_map.keys())
    aa_pattern = '|'.join(aa_codes)
    
    # Try three-letter format first
    pattern = rf'^({aa_pattern})([0-9]+)({aa_pattern})$'
    match = re.match(pattern, clean_seq, re.IGNORECASE)
    if match:
        ref_aa_3 = match.group(1).capitalize()
        position = int(match.group(2))
        alt_aa_3 = match.group(3).capitalize()
        
        ref_aa = aa_map.get(ref_aa_3, ref_aa_3[0])  # Fallback to first letter
        alt_aa = aa_map.get(alt_aa_3, alt_aa_3[0])
        
        return ref_aa, position, alt_aa
    
    # Try single-letter format
    pattern = r'^([A-Z*X])([0-9]+)([A-Z*X])$'
    match = re.match(pattern, clean_seq, re.IGNORECASE)
    if match:
        ref_aa = match.group(1).upper()
        position = int(match.group(2))
        alt_aa = match.group(3).upper()
        return ref_aa, position, alt_aa
    
    raise ValueError(f"Invalid AA change format: {variant_seq}")

def detect_variant_format(counts: Dict[int, Dict[int, pd.DataFrame]]) -> str:
    """
    Auto-detect the format of variant sequences in count data.
    
    Parameters
    ----------
    counts : dict
        Nested dict of count DataFrames
        
    Returns
    -------
    str
        'aa_changes' if pre-annotated AA changes detected, 'aa' for full sequences, 'dna' for DNA
    """
    # Get a sample of variant sequences
    sample_variants = []
    for rep_dict in counts.values():
        for df in rep_dict.values():
            sample_variants.extend(df['variant_seq'].head(20).tolist())
            if len(sample_variants) >= 50:
                break
        if len(sample_variants) >= 50:
            break
    
    # Check if most variants match AA change format
    aa_change_count = sum(1 for seq in sample_variants if is_aa_change_format(seq))
    
    if aa_change_count / len(sample_variants) > 0.8:  # 80% threshold
        return 'aa_changes'
    
    # Check average length to distinguish between full AA sequences and DNA
    avg_length = sum(len(seq) for seq in sample_variants) / len(sample_variants)
    
    if avg_length > 100:  # Likely DNA sequences
        return 'dna'
    else:
        return 'aa'  # Full AA sequences

@dataclass
class ExperimentConfig:
    """
    Dataclass for experiment configuration and data loading.
    """
    experiment_name: str
    experiment_setup_file: str
    wt_seq: str
    analysis_type: str  # 'snv', 'codon', or 'aa'
    output_dir: Optional[str] = None
    other_params: Optional[Dict[str, Any]] = None
    counts: Optional[Dict[int, Dict[int, pd.DataFrame]]] = None
    mfi: Optional[Dict[int, Dict[int, float]]] = None
    total_reads: Optional[Dict[int, Dict[int, int]]] = None
    cell_prop: Optional[Dict[int, Dict[int, float]]] = None
    
    # Auto-detected properties (set during loading)
    variant_type: Optional[str] = None  # Auto-detected: 'dna' or 'aa'
    min_pos: Optional[int] = None  # Auto-detected from data
    max_pos: Optional[int] = None  # Auto-detected from data
    position_type: str = 'aa'  # Controls heatmap x-axis ('aa' or 'dna')
    
    # Analysis parameters with defaults
    bins_required: int = 1
    reps_required: int = 1  
    avg_method: str = 'rep-weighted'
    max_cv: Optional[float] = None
    minread_threshold: int = 0
    barcoded: bool = False
    aa_pre_annotated: bool = False
    transcript_id: Optional[str] = None
    mutagenesis_variants: Optional[list] = None
    position_offset: int = 0  # Offset for position numbering (e.g., if data positions start from 1 but gene positions start from 51)
    biophysical_prop: bool = False  # Whether to show biophysical properties panel in heatmaps
    
    @property
    def num_aa(self) -> int:
        """Calculate number of amino acids from detected position range."""
        if self.min_pos is None or self.max_pos is None:
            self._detect_position_range()
        return self.max_pos - self.min_pos + 1
    
    @property 
    def num_positions(self) -> int:
        """Calculate number of positions based on variant_type."""
        if self.variant_type == 'dna':
            # For DNA positions, use the full DNA sequence length
            return len(self.wt_seq)
        else:
            # For AA positions, use the detected range
            if self.min_pos is None or self.max_pos is None:
                self._detect_position_range()
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
        
        # Validate analysis_type
        valid_analysis_types = {'snv', 'codon', 'aa'}
        if 'analysis_type' not in data:
            raise ValueError(f"Missing required 'analysis_type' parameter. Must be one of: {valid_analysis_types}")
        if data['analysis_type'] not in valid_analysis_types:
            raise ValueError(f"Invalid analysis_type '{data['analysis_type']}'. Must be one of: {valid_analysis_types}")

        # Required fields
        for field in ['experiment_name', 'experiment_setup_file', 'wt_seq', 'analysis_type']:
            if field in data:
                args[field] = data[field]
        
        # Optional fields (only add if present in JSON to preserve dataclass defaults)
        optional_fields = [
            'output_dir', 'bins_required', 'reps_required', 'avg_method',
            'minread_threshold', 'barcoded', 'max_cv', 'transcript_id',
            'mutagenesis_variants', 'position_offset', 'biophysical_prop',
            'position_type'
        ]

        for field in optional_fields:
            if field in data:
                args[field] = data[field]
        
        # Other parameters
        handled_keys = {'experiment_name', 'experiment_setup_file', 'wt_seq', 'analysis_type',
                       'output_dir', 'bins_required', 'reps_required', 'avg_method', 'minread_threshold', 'barcoded', 'max_cv', 'transcript_id', 'mutagenesis_variants', 'position_offset', 'biophysical_prop', 'position_type'}
        other_params = {k: v for k, v in data.items() if k not in handled_keys}
        if other_params:
            args['other_params'] = other_params
        
        # Create config instance
        config = ExperimentConfig(**args)
        
        # Auto-detect variant_type from count files
        from sortscore.analysis.variant_detection import detect_variant_type_from_experiment
        try:
            config.variant_type = detect_variant_type_from_experiment(config.experiment_setup_file)
            logging.info(f"Auto-detected variant_type: '{config.variant_type}'")
        except Exception as e:
            raise ValueError(f"Failed to auto-detect variant_type from count files: {e}")

        # Set default position_type when not provided: SNV analyses default to DNA coordinates
        position_type = data.get('position_type')
        if position_type is None:
            config.position_type = 'dna' if config.analysis_type == 'snv' else 'aa'
        elif position_type not in {'aa', 'dna'}:
            raise ValueError("position_type must be 'aa' or 'dna'")
        else:
            config.position_type = position_type
        
        # Load counts to enable position detection
        config.load_counts()
        
        # Validate analysis_type compatibility with detected variant_type
        if config.analysis_type == 'aa' and config.variant_type == 'dna':
            logging.info("Analysis type 'aa' with DNA variants - will aggregate synonymous variants")
        elif config.analysis_type in ['snv', 'codon'] and config.variant_type == 'aa':
            raise ValueError(f"Analysis type '{config.analysis_type}' requires DNA variants, but detected AA variants")
        elif config.analysis_type == 'aa' and config.variant_type == 'aa':
            logging.info("Analysis type 'aa' with AA variants - direct processing")
        
        return config

    def load_counts(self) -> None:
        """
        Load count DataFrames for all replicates and bins as specified in the experiment setup file.

        This populates the `counts` attribute as a nested dictionary: counts[rep][bin] = DataFrame.
        Also populates the `mfi` attribute as a nested dictionary: mfi[rep][bin] = float.

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
        mfi = {}
        total_reads = {}
        cell_prop = {}
        
        for _, row in setup_df.iterrows():
            rep = int(row['Replicate'])
            bin_ = str(row['Bin']).strip()
            count_file = str(row['Path']).strip()
            if 'Read Counts (CSV)' in row:
                count_file = str(row['Read Counts (CSV)']).strip()
                
            mfi_val = float(row['MFI'])
            
            # Load total reads if available (for sample read depth normalization)
            if 'Read Count' in row:
                total_read_count = int(row['Read Count'])
                total_reads.setdefault(rep, {})[bin_] = total_read_count
            
            # Load cell proportions if available (for cell distribution normalization)
            if 'Proportion of Cells' in row:
                cell_proportion = float(row['Proportion of Cells'])
                cell_prop.setdefault(rep, {})[bin_] = cell_proportion
            
            try:
                # Check if file is parquet format
                if count_file.endswith('.parquet'):
                    df = self._load_parquet_as_counts(count_file)
                else:
                    # Auto-detect header by checking if second column can be converted to numeric
                    header_row = self._detect_header_row(count_file)
                    df = pd.read_csv(count_file, sep=None, engine='python', header=header_row)
                    
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
            mfi.setdefault(rep, {})[bin_] = mfi_val
            
        self.counts = counts
        self.mfi = mfi
        
        # Auto-detect variant format and convert pre-annotated AA changes if needed
        if counts:
            detected_format = detect_variant_format(counts)
            if detected_format == 'aa_changes':
                logging.info("Detected pre-annotated AA changes format. Converting to full annotations...")
                self._convert_aa_changes_to_annotations()
        
        # Set total_reads if we loaded any
        if total_reads:
            self.total_reads = total_reads
            
        # Set cell_prop if we loaded any
        if cell_prop:
            self.cell_prop = cell_prop
        
        # Auto-detect position range from loaded data
        self._detect_position_range()
    
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
    
    def get_merged_counts(self) -> pd.DataFrame:
        """
        Create a merged DataFrame with all variants and their counts across replicates/bins.
        
        This preserves all columns from the original count DataFrames (like aa_seq_diff)
        and creates a wide-format DataFrame with count columns for each rep/bin.
        
        Returns
        -------
        pd.DataFrame
            Merged DataFrame with variant_seq and additional columns preserved,
            plus count columns for each replicate/bin combination.
        """
        if not self.counts:
            raise ValueError("No counts loaded. Call load_counts() first.")
        
        # Get all unique variants and collect additional columns from first available DataFrame
        all_variants = set()
        additional_columns = {}
        first_df = None
        
        for rep in self.counts:
            for bin_ in self.counts[rep]:
                all_variants.update(self.counts[rep][bin_]['variant_seq'])
                if first_df is None:
                    first_df = self.counts[rep][bin_]
                    # Get any additional columns beyond variant_seq and count
                    for col in first_df.columns:
                        if col not in ['variant_seq', 'count']:
                            additional_columns[col] = first_df.set_index('variant_seq')[col].to_dict()
        
        all_variants = sorted(all_variants)
        df = pd.DataFrame({'variant_seq': all_variants})
        
        # Add additional columns (like aa_seq_diff) if they exist
        for col_name, col_data in additional_columns.items():
            df[col_name] = df['variant_seq'].map(col_data).fillna('')
        
        # Add count columns for each rep/bin
        for rep in self.counts:
            for bin_ in self.counts[rep]:
                col = f'count.r{rep}b{bin_}'
                d = self.counts[rep][bin_].set_index('variant_seq')['count']
                df[col] = df['variant_seq'].map(d).fillna(0)
        
        return df

    def _load_parquet_as_counts(self, parquet_file: str) -> pd.DataFrame:
        """
        Load parquet file and convert to variant counts format.
        
        Expected parquet format:
        - Contains columns: var_ref, var_pos, var_alt, read_count
        - Creates variant identifiers like "A36T" from ref/pos/alt
        - Aggregates read counts by variant
        
        Parameters
        ----------
        parquet_file : str
            Path to parquet file
            
        Returns
        -------
        pd.DataFrame
            DataFrame with columns ['variant_seq', 'count']
        """
        # Read parquet file
        df = pd.read_parquet(parquet_file)
        
        # Filter to only mapped variants (non-null var_pos)
        mapped_df = df[df['var_pos'].notna()].copy()
        
        # Create variant identifier (e.g., "A36T", "G252M")
        mapped_df['variant'] = (mapped_df['var_ref'].astype(str) + 
                               mapped_df['var_pos'].astype(int).astype(str) + 
                               mapped_df['var_alt'].astype(str))
        
        # Aggregate read counts by variant
        variant_counts = mapped_df.groupby('variant')['read_count'].sum().reset_index()
        
        # Rename columns for sortscore compatibility
        variant_counts.columns = ['variant_seq', 'count']
        
        # Sort by count (descending)
        variant_counts = variant_counts.sort_values('count', ascending=False)
        
        logging.info(f"Loaded parquet file {parquet_file}: {len(variant_counts)} variants, {variant_counts['count'].sum():,} total reads")
        
        return variant_counts

    def _detect_header_row(self, count_file: str) -> int:
        """
        Auto-detect the header row by finding the first line where the second column 
        can be successfully converted to numeric.
        
        Parameters
        ----------
        count_file : str
            Path to the count file
            
        Returns
        -------
        int
            Row number to use as header (0-based), or None if no header detected
        """
        try:
            import gzip
            # Handle both regular and gzipped files
            if count_file.endswith('.gz'):
                open_func = lambda f: gzip.open(f, 'rt')
            else:
                open_func = lambda f: open(f, 'r')
                
            # Read first few lines to check for headers
            with open_func(count_file) as f:
                lines = []
                for _ in range(5):  # Read up to 5 lines
                    line = f.readline()
                    if not line:
                        break
                    lines.append(line.strip())
            
            # Test each line to see if second column is numeric
            for line_idx, line in enumerate(lines):
                if not line:
                    continue
                    
                # Try to parse the line with common separators
                for sep in ['\t', ',', ' ']:
                    try:
                        parts = line.split(sep)
                        if len(parts) >= 2:
                            # Try to convert second column to numeric
                            pd.to_numeric(parts[1].strip())
                            # If successful, this line contains data, so header is the line before
                            return line_idx - 1 if line_idx > 0 else None
                    except (ValueError, pd.errors.ParserError):
                        continue
            
            # If no numeric second column found, assume no header
            return None
            
        except Exception:
            # If file reading fails, default to header=0
            return 0

    def _convert_aa_changes_to_annotations(self) -> None:
        """
        Convert pre-annotated AA changes to the format expected by existing annotation functions.
        
        This creates an aa_seq_diff column in the format 'ref.position.alt' for differences
        from the wild-type sequence, and uses the existing annotation function to determine annotation types.
        """
        from sortscore.sequence_parsing import translate_dna
        
        # Get the wild-type AA sequence for annotation
        wt_aa_seq = translate_dna(self.wt_seq)
        
        for rep_dict in self.counts.values():
            for df in rep_dict.values():
                aa_seq_diffs = []
                
                for variant_seq in df['variant_seq']:
                    try:
                        ref_aa, position, alt_aa = parse_aa_change(variant_seq)
                        
                        # Verify the reference AA matches the WT sequence
                        if position <= len(wt_aa_seq) and wt_aa_seq[position-1] == ref_aa:
                            # Only add difference if alt_aa is different from ref_aa
                            if ref_aa != alt_aa:
                                aa_seq_diff = f"{ref_aa}.{position}.{alt_aa}"
                                aa_seq_diffs.append(aa_seq_diff)
                            else:
                                # No difference (e.g., M1M)
                                aa_seq_diffs.append('')
                        else:
                            # Reference doesn't match WT, still record the difference
                            aa_seq_diff = f"{ref_aa}.{position}.{alt_aa}"
                            aa_seq_diffs.append(aa_seq_diff)
                        
                    except ValueError:
                        # If parsing fails, treat as empty difference
                        aa_seq_diffs.append('')
                
                # Add the aa_seq_diff column
                df['aa_seq_diff'] = aa_seq_diffs
                
                # Add annotation categories
                from sortscore.analysis.annotation import add_variant_categories
                df = add_variant_categories(df)
    
    def _detect_position_range(self) -> None:
        """Auto-detect min_pos and max_pos from loaded variant data using sequence comparison."""
        if not self.counts:
            logging.warning("No counts loaded, cannot detect position range")
            return
            
        min_pos = float('inf')
        max_pos = float('-inf')
        
        try:
            ref_seq = get_reference_sequence(self.wt_seq, self.variant_type)
        except ValueError as e:
            logging.error(f"Cannot get reference wild-type sequence: {e}")
            return


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
