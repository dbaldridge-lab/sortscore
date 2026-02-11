"""
Automatic variant type detection from count files.

This module provides functions to automatically detect the format of variant
sequences in count files. It checks for common variant notations using HGVS 
parsing and pattern matching. For full sequences, pattern matching is used to
distinguish between sequence types. This module can detect mutagenesis types such
as amino acid (AA) substitutions, codon changes, and single nucleotide variants (SNV).
"""
import pandas as pd
import re
from typing import List, Tuple, Literal, Optional
from pathlib import Path

from sortscore.utils.experiment_setup import load_experiment_setup


def detect_sequence_format(sequences: List[str]) -> Literal['dna', 'aa']:
    """
    Classify a set of sequences as containing DNA or AA input sequences
    using HGVS parsing and pattern matching.
    
    Analyzes a sample of sequences from an input file to determine if they are:
    - DNA sequences (full sequences or nucleotide changes)
    - Amino acid sequences (full sequences or AA changes)
    
    Parameters
    ----------
    sequences : List[str]
        Sample of sequence strings to analyze
        
    Returns
    -------
    Literal['dna', 'aa']
        Detected sequence format
        
    Examples
    --------
    # TODO: add DNA variant and ID examples T123G, c.123A>G, g.123456A>G
    # TODO: could someone provide multiple DNA variants in the input column
    >>> sequences = ['ATGCGT', 'ATGCGA', 'ATGCGG']
    >>> detect_sequence_format(sequences)
    'dna'
    
    # TODO: add AA sequence and ID examples
    >>> sequences = ['M1V', 'R98C', 'P171*']
    >>> detect_sequence_format(sequences)
    'aa'
    
    Raises
    ------
    ValueError
        If mixed formats are detected or format cannot be determined
    """
    if not sequences:
        raise ValueError("No sequences provided for format detection")
    
    # Remove empty/null sequences
    valid_sequences = [seq for seq in sequences if seq and pd.notna(seq) and str(seq).strip()]
    if not valid_sequences:
        raise ValueError("No valid sequences provided for format detection")

    dna_count = 0
    aa_count = 0
    
    for seq in valid_sequences:
        seq = str(seq).strip()
        if not seq:
            raise ValueError("Empty sequence encountered during format detection")
            
        format_type = _classify_single_sequence(seq)
        if format_type == 'dna':
            dna_count += 1
        elif format_type == 'aa':
            aa_count += 1
    
    # Determine overall format
    total_classified = dna_count + aa_count
    if total_classified == 0:
        raise ValueError("All input sequences are unrecognized. Could not determine sequence format.")
    
    dna_fraction = dna_count / total_classified
    aa_fraction = aa_count / total_classified
    
    # Require consensus to avoid mixed format issues
    if dna_fraction == 1.0:
        return 'dna'
    elif aa_fraction == 1.0:
        return 'aa'
    elif dna_count > 0 and aa_count > 0:
        raise ValueError("Mixed sequence formats detected (both DNA and AA present). Ensure all sequences use a consistent format.")
    else:
        raise ValueError("Could not determine sequence format from input sequences.")

def _parse_hgvs_variant(seq: str) -> Optional[Literal['dna', 'aa']]:
    """
    Try to parse sequence as HGVS notation and return variant type.
    Parameters
    ----------
    seq : str
        Sequence to parse as HGVS
    Returns
    -------
    Optional[Literal['dna', 'aa']]
        'aa' for protein variants (p.), 'dna' for DNA variants (c./g./n.), None if not HGVS
    """
    try:
        from mavehgvs import Variant

        variant = Variant(seq)

        # Check variant type based on HGVS prefix
        if variant.prefix == 'p':
            return 'aa'  # Protein variant
        elif variant.prefix in ['c', 'g', 'n']:
            return 'dna'  # DNA variant
        else:
            return None

    except ImportError:
        # mavehgvs library not available, fall back to pattern matching
        return None
    except Exception:
        # Not valid HGVS notation or parsing error
        return None


def _classify_single_sequence(seq: str) -> Optional[Literal['dna', 'aa']]:
    """
    Classify a single sequence as DNA, AA, or unknown.

    Uses HGVS parsing first, then falls back to pattern matching.

    Parameters
    ----------
    seq : str
        Sequence string to classify

    Returns
    -------
    Literal['dna', 'aa', 'unknown']
        Classification of the sequence
    """
    seq = seq.strip()

    # Try HGVS parsing (case-sensitive)
    hgvs_type = _parse_hgvs_variant(seq)
    if hgvs_type:
        return hgvs_type

    # For pattern matching, use uppercase
    seq = seq.upper()

    # Fall back to pattern matching for common non-HGVS formats
    
    # Amino acid change patterns (single letter)
    aa_single_patterns = [
        r'^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY\*]$',  # M1V, R98*
        r'^[ACDEFGHIKLMNPQRSTVWY]\.\d+\.[ACDEFGHIKLMNPQRSTVWY\*]$',  # M.1.V, R.98.*
        r'^[ACDEFGHIKLMNPQRSTVWY]-\d+-[ACDEFGHIKLMNPQRSTVWY\*]$',  # M-1-V, R-98-*
        r'^[ACDEFGHIKLMNPQRSTVWY]_\d+_[ACDEFGHIKLMNPQRSTVWY\*]$',  # M_1_V, R_98_*
    ]
    
    for pattern in aa_single_patterns:
        if re.match(pattern, seq):
            return 'aa'
    
    # Amino acid change patterns (three letter)
    aa_three_letter = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'TER'
    ]
    
    # Three letter AA change patterns
    three_letter_pattern = r'^([A-Z]{3})(\d+)([A-Z]{3})$'
    match = re.match(three_letter_pattern, seq)
    if match:
        from_aa, pos, to_aa = match.groups()
        if from_aa in aa_three_letter and to_aa in aa_three_letter:
            return 'aa'
    
    # DNA change patterns (nucleotide substitutions)
    dna_change_patterns = [
        r'^[ATCG]\d+[ATCG]$',  # A123T
        r'^[ATCG]\.\d+\.[ATCG]$',  # A.123.T
        r'^[ATCG]-\d+-[ATCG]$',  # A-123-T
        r'^[ATCG]_\d+_[ATCG]$',  # A_123_T
    ]
    
    for pattern in dna_change_patterns:
        if re.match(pattern, seq):
            return 'dna'
    
    # Full sequence detection
    dna_bases = set('ATCG')
    aa_codes = set('ACDEFGHIKLMNPQRSTVWY*')
    
    seq_chars = set(seq)
    
    # If sequence contains only DNA bases
    if seq_chars.issubset(dna_bases) and len(seq) > 1:
        return 'dna'
    
    # If sequence contains only amino acid codes
    if seq_chars.issubset(aa_codes) and len(seq) > 1:
        return 'aa'
    
    # Check for mixed bases that suggest DNA (includes N, ambiguous bases)
    extended_dna_bases = set('ATCGRYSWKMBDHVN')
    if seq_chars.issubset(extended_dna_bases) and len(seq) > 1:
        return 'dna'
    
    return None

# TODO: detect codon and snv variant logic
# TODO: Pass filename to detect_sequence_format
def detect_variant_type_from_counts(counts) -> str:
    """
    Detect variant type by examining loaded count DataFrames.

    Parameters
    ----------
    counts : dict
        Nested dict: counts[rep][bin] = DataFrame with a 'variant_seq' column.

    Returns
    -------
    str
        Detected variant type ('dna' or 'aa').

    Raises
    ------
    ValueError
        If no sequences are found, mixed formats are detected, or sequences are unrecognized.
    """
    if not counts:
        raise ValueError("No counts provided for variant type detection")

    dna_count = 0
    aa_count = 0
    unknown_examples: List[str] = []

    for rep_dict in counts.values():
        for df in rep_dict.values():
            if df is None or df.empty or 'variant_seq' not in df.columns:
                continue
            for seq in df['variant_seq'].dropna().astype(str):
                seq = seq.strip()
                if not seq:
                    continue
                format_type = _classify_single_sequence(seq)
                if format_type == 'dna':
                    dna_count += 1
                elif format_type == 'aa':
                    aa_count += 1
                else:
                    if len(unknown_examples) < 5:
                        unknown_examples.append(seq)
                if dna_count > 0 and aa_count > 0:
                    raise ValueError(
                        "Mixed sequence formats detected (both DNA and AA present). "
                        "Ensure all sequences use a consistent format."
                    )

    if dna_count == 0 and aa_count == 0:
        raise ValueError("No valid sequences found in loaded counts for variant type detection")

    if unknown_examples:
        raise ValueError(
            "Unrecognized sequence format(s) detected in count files. "
            f"Examples: {unknown_examples}"
        )

    return 'dna' if dna_count > 0 else 'aa'


def detect_variant_type_from_experiment(experiment_setup_file: str) -> str:
    """
    Detect variant type by examining count files from experiment setup. 
    
    Parameters
    ----------
    experiment_setup_file : str
        Path to experiment setup CSV file
        
    Returns
    -------
    str
        Detected variant type ('aa', 'snv', 'codon')
        
    Raises
    ------
    ValueError
        If variant type cannot be determined or is mixed/ambiguous
    """
    setup_df, setup_cols = load_experiment_setup(experiment_setup_file)
    
    # Sample sequences from all count files to ensure consistency
    all_sequences = []
    files_sampled = 0
    setup_dir = Path(str(experiment_setup_file).strip()).expanduser().resolve().parent
    
    for file_path in setup_df[setup_cols.count_file].dropna():
        try:
            file_path = setup_dir / Path(str(file_path).strip()).expanduser()
            file_path = file_path.resolve()
            
            # Load read counts from file
            count_df = pd.read_csv(file_path, sep=None, engine='python')
            
            # Get sequence column (assumes first column is sequences)
            seq_col = count_df.columns[0]
            sequences = count_df[seq_col].dropna().astype(str).tolist()
            
            # Take sample from the file
            all_sequences.extend(sequences[:100])
            files_sampled += 1
            
        except Exception as e:
            raise RuntimeError(f"Failed to process count file '{file_path}': {e}") from e
    
    if not all_sequences:
        raise ValueError(
            f"Could not read any sequences from count files listed in setup file: {experiment_setup_file}"
        )

    detected_format = detect_sequence_format(all_sequences)

    return detected_format
