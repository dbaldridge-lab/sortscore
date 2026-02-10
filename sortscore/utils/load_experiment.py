"""
Configuration and constants for Sort-seq variant analysis.

This module provides the ExperimentConfig dataclass for experiment configuration and data loading.
It includes automatic detection and conversion of pre-annotated amino acid changes.

Supported Variant Formats
-------------------------
1. **Full Sequences**: DNA or protein sequences to be compared against the wild-type sequence.
   - DNA: Full nucleotide sequences for variant_type="codon" or "snv"
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
from pathlib import Path
import json
import logging
import pandas as pd
import re
from dataclasses import dataclass
from typing import Dict, Any, Optional

from sortscore.utils.sequence_parsing import get_reference_sequence
from sortscore.utils.experiment_setup import load_experiment_setup

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
    output_dir: Optional[str] = None
    other_params: Optional[Dict[str, Any]] = None
    counts: Optional[Dict[int, Dict[int, pd.DataFrame]]] = None
    mfi: Optional[Dict[int, Dict[int, float]]] = None
    total_reads: Optional[Dict[int, Dict[int, int]]] = None
    cell_prop: Optional[Dict[int, Dict[int, float]]] = None
    position_type: str = 'aa'  # Controls heatmap x-axis ('aa' or 'dna')

    # Auto-detected properties (set during loading)
    variant_type: Optional[str] = None  # Auto-detected: 'aa', 'codon', 'snv'"
    min_pos: Optional[int] = None  # Auto-detected from data
    max_pos: Optional[int] = None  # Auto-detected from data

    # Analysis parameters with defaults
    bins_required: int = 1
    reps_required: int = 1  
    avg_method: str = 'rep-weighted'
    max_cv: Optional[float] = None
    minread_threshold: int = 0
    aa_pre_annotated: bool = False
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
        if self.variant_type == 'snv':
            # For DNA positions, use the full DNA sequence length
            return len(self.wt_seq)
        if self.variant_type == 'codon':
            return len(self.wt_seq) // 3
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
        json_path_obj = Path(json_path).expanduser().resolve()
        with open(json_path_obj, 'r') as f:
            data = json.load(f)
        return ExperimentConfig.from_dict(data, config_file_dir=json_path_obj.parent)

    @staticmethod
    def from_dict(data: Dict[str, Any], config_file_dir: Optional[Path] = None) -> 'ExperimentConfig':
        """
        Load experiment configuration from an in-memory mapping.

        Implements CLI-based configuration (merged from CLI args + optional JSON config).

        Parameters
        ----------
        data : Dict[str, Any]
            Mapping of configuration keys to values (e.g. parsed JSON).
        config_file_dir : Optional[pathlib.Path]
            Base directory for resolving relative paths found in `data`.
            If None, relative paths are resolved from the current working directory.
        """
        # Build kwargs only for fields that exist in data
        args: Dict[str, Any] = {}
        base_dir = config_file_dir or Path.cwd()

        # Required fields
        required_fields = ['experiment_name', 'experiment_setup_file', 'wt_seq']
        missing_required = [field for field in required_fields if field not in data or data[field] in (None, "")]
        if missing_required:
            raise ValueError(
                "Missing required configuration value(s): "
                + ", ".join(missing_required)
                + ". Provide them via CLI arguments or an optional JSON config (-c)."
            )

        for field in required_fields:
            value = data[field]
            if field == 'experiment_setup_file':
                value = str((base_dir / Path(str(value)).expanduser()).resolve())
            args[field] = value

        # Optional fields (only add if present in JSON to preserve dataclass defaults)
        optional_fields = [
            'output_dir', 'bins_required', 'reps_required', 'avg_method',
            'minread_threshold','max_cv',
            'mutagenesis_variants', 'position_offset', 'biophysical_prop',
            'position_type', 'min_pos', 'max_pos'
        ]

        for field, value in data.items():
            if field in optional_fields:
                if field == 'output_dir' and value is not None:
                    value = str((base_dir / Path(str(value)).expanduser()).resolve())
                args[field] = value
        
        # Other parameters
        handled_keys = {'experiment_name', 'experiment_setup_file', 'wt_seq', 
                       'output_dir', 'bins_required', 'reps_required', 'avg_method', 
                       'minread_threshold', 'max_cv', 'mutagenesis_variants', 
                       'position_offset', 'biophysical_prop', 'position_type'}
        other_params = {k: v for k, v in data.items() if k not in handled_keys}
        if other_params:
            args['other_params'] = other_params
        
        # Create config instance
        config = ExperimentConfig(**args)

        # Validate experiment setup CSV early for clearer errors
        try:
            load_experiment_setup(config.experiment_setup_file)
        except Exception as e:
            raise ValueError(f"Invalid experiment setup CSV '{config.experiment_setup_file}': {e}") from e

        # Load counts once; variant type detection uses loaded sequences (no extra file reads).
        config.load_counts(detect_position_range=False)

        from sortscore.utils.variant_detection import detect_variant_type_from_counts
        try:
            config.variant_type = detect_variant_type_from_counts(config.counts)
            logging.info(f"Auto-detected variant_type: '{config.variant_type}'")
        except Exception as e:
            raise ValueError(f"Failed to auto-detect variant_type from loaded counts: {e}") from e

        # Set default position_type when not provided: SNV analyses default to DNA coordinates
        position_type = data.get('position_type')
        if position_type is None:
            config.position_type = 'dna' if config.variant_type == 'snv' else 'aa'
        elif position_type not in {'aa', 'dna'}:
            raise ValueError("position_type must be 'aa' or 'dna'")
        else:
            config.position_type = position_type

        # Now that variant_type/position_type are set, detect position range from loaded data.
        config._detect_position_range()
        
        return config

    def load_counts(self, detect_position_range: bool = True) -> None:
        """
        Load count DataFrames for all replicates and bins as specified in the experiment setup file.

        This populates the `counts` attribute as a nested dictionary: counts[rep][bin] = DataFrame.
        Also populates the `mfi` attribute as a nested dictionary: mfi[rep][bin] = float.

        Returns
        -------
        None
        """
        try:
            setup_df, setup_cols = load_experiment_setup(self.experiment_setup_file)
        except Exception as e:
            raise RuntimeError(f"Failed to load experiment setup file: {e}") from e
        counts = {}
        mfi = {}
        total_reads = {}
        cell_prop = {}
        config_file_dir = Path(self.experiment_setup_file).expanduser().resolve().parent
        
        for _, row in setup_df.iterrows():
            rep = int(row[setup_cols.replicate])
            bin_ = str(row[setup_cols.bin]).strip()
            count_file = str(row[setup_cols.count_file]).strip()
            count_file = config_file_dir / Path(count_file).expanduser()
            count_file = count_file.resolve()
            count_file = str(count_file)
            
            mfi_val = float(row[setup_cols.mfi])
            
            # Load total reads if available (for sample read depth normalization)
            if setup_cols.read_count is not None:
                read_count_val = row[setup_cols.read_count]
                if pd.notna(read_count_val):
                    total_read_count = int(read_count_val)
                    total_reads.setdefault(rep, {})[bin_] = total_read_count
            
            # Load cell proportions if available (for cell distribution normalization)
            if setup_cols.proportion_of_cells is not None:
                cell_prop_val = row[setup_cols.proportion_of_cells]
                if pd.notna(cell_prop_val):
                    cell_proportion = float(cell_prop_val)
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
            
        # Fail fast if nothing usable was loaded (prevents later redundant errors in variant detection).
        has_any_sequences = False
        for rep_dict in counts.values():
            for df in rep_dict.values():
                if df is None or df.empty or 'variant_seq' not in df.columns:
                    continue
                if df['variant_seq'].dropna().astype(str).str.strip().ne('').any():
                    has_any_sequences = True
                    break
            if has_any_sequences:
                break
        if not has_any_sequences:
            raise ValueError(
                f"Could not read any sequences from count files listed in setup file: {self.experiment_setup_file}"
            )

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
        if detect_position_range:
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
        from sortscore.utils.sequence_parsing import translate_dna
        
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
        # If both values are already set (e.g., provided in config), keep them
        if self.min_pos is not None and self.max_pos is not None:
            return

        if not self.counts:
            logging.warning("No counts loaded, cannot detect position range")
            # Fallback to sequence length so downstream code has usable values
            try:
                ref_seq = get_reference_sequence(self.wt_seq, self.variant_type)
                seq_len = len(ref_seq)
            except Exception:
                seq_len = len(self.wt_seq)
            self.min_pos = self.min_pos or 1
            self.max_pos = self.max_pos or seq_len
            return
            
        min_pos = float('inf')
        max_pos = float('-inf')
        
        try:
            ref_seq = get_reference_sequence(self.wt_seq, self.variant_type)
        except ValueError as e:
            logging.error(f"Cannot get reference wild-type sequence: {e}")
            ref_seq = self.wt_seq

        # If no variant annotations are available, fall back to full sequence length
        if self.min_pos is None:
            self.min_pos = 1
        if self.max_pos is None:
            # TODO: will this work with variant IDs that are not full sequences? May need to edit ref_seq
            if self.variant_type in {'codon', 'snv', 'dna'} and self.position_type == 'dna':
                self.max_pos = len(ref_seq)
            else:
                aa_len = len(get_reference_sequence(ref_seq, 'aa')) if ref_seq else 0
                self.max_pos = aa_len
        # If values were partially provided (one missing), fill the other based on sequence length
        if self.min_pos is None:
            self.min_pos = 1
        if self.max_pos is None:
            self.max_pos = len(ref_seq)

    def annotate_variants(self, df, wt_ref_seq: str, wt_aa_seq: str, variant_type: str) -> None:
        """
        Annotate a variant count DataFrame with sequence difference from WT.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to annotate (will be modified in place).
        wt_ref_seq : str
            Wild-type DNA sequence.
        wt_aa_seq : str
            Wild-type amino acid sequence.
        variant_type : str
            'aa' for amino acid variants, 'codon' for multiple nucleotide changes in single frame, 'snv' for single nucleotide variants.
        """
        from sortscore.utils.sequence_parsing import compare_to_reference, translate_dna
        
        # AA annotation
        if variant_type == 'aa':
            df['aa_seq'] = df['variant_seq']  # AA sequences provided directly
        # TODO: will this work with variant IDs that are not full sequences?
        elif variant_type in {'codon', 'snv', 'dna'}:
            df['aa_seq'] = df['variant_seq'].apply(translate_dna)  # Translate DNA to AA
            # DNA annotation
            df['dna_seq_diff'] = df['variant_seq'].apply(lambda x: compare_to_reference(wt_ref_seq, x))
            df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
            df['dna_diff_count'] = df['dna_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
            
        df['aa_seq_diff'] = df['aa_seq'].apply(lambda x: compare_to_reference(wt_aa_seq, x))
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        df['aa_diff_count'] = df['aa_seq_diff'].apply(lambda x: 0 if not x else x.count(',') + 1)
