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
import logging

from sortscore.analysis.experiment_setup import load_experiment_setup


def detect_sequence_format(sequences: List[str]) -> Literal['dna', 'aa', 'mixed', 'unknown']:
    """
    Destinguish between DNA and AA input sequences using HGVS parsing and pattern matching.
    
    Analyzes a sample of sequences from an input file to determine if they are:
    - DNA sequences (full sequences or nucleotide changes)
    - Amino acid sequences (full sequences or AA changes) 
    - Mixed formats (error condition)
    - Unknown format (error condition)
    
    Parameters
    ----------
    sequences : List[str]
        Sample of sequence strings to analyze
        
    Returns
    -------
    Literal['dna', 'aa', 'mixed', 'unknown']
        Detected sequence format
        
    Examples
    --------
    # TODO: add DNA variant examples T123G, c.123A>G, g.123456A>G
    # TODO: could someone provide multiple DNA variants in the input column
    >>> sequences = ['ATGCGT', 'ATGCGA', 'ATGCGG']
    >>> detect_sequence_format(sequences)
    'dna'
    
    # TODO: add AA sequence examples
    >>> sequences = ['M1V', 'R98C', 'P171*']
    >>> detect_sequence_format(sequences)
    'aa'
    """
    if not sequences:
        return 'unknown'
    
    # Remove empty/null sequences
    valid_sequences = [seq for seq in sequences if seq and pd.notna(seq) and str(seq).strip()]
    if not valid_sequences:
        return 'unknown'
    
    # Take a sample for analysis (max 100 sequences)
    sample = valid_sequences[:100]
    
    dna_count = 0
    aa_count = 0
    unknown_count = 0
    
    for seq in sample:
        seq = str(seq).strip()
        if not seq:
            continue
            
        format_type = _classify_single_sequence(seq)
        if format_type == 'dna':
            dna_count += 1
        elif format_type == 'aa':
            aa_count += 1
        else:
            unknown_count += 1
    
    # Determine overall format
    total_classified = dna_count + aa_count + unknown_count
    if total_classified == 0:
        return 'unknown'
    
    dna_fraction = dna_count / total_classified
    aa_fraction = aa_count / total_classified
    unknown_fraction = unknown_count / total_classified
    
    logging.info(
        f"Sequence format fractions: DNA={dna_fraction:.2f}, AA={aa_fraction:.2f}, Unknown={unknown_fraction:.2f}"
    )
    
    # Require consensus to avoid mixed format issues
    if dna_fraction == 1.0:
        return 'dna'
    elif aa_fraction == 1.0:
        return 'aa'
    else:
        return 'mixed'


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


def _classify_single_sequence(seq: str) -> Literal['dna', 'aa', 'unknown']:
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
    
    return 'unknown'

# TODO: detect codon and snv variant logic
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
            # Handle paths relative to setup file
            file_path = setup_dir / Path(str(file_path).strip()).expanduser()
            file_path = file_path.resolve()
            
            # Read count file
            count_df = pd.read_csv(file_path, sep=None, engine='python')  # Auto-detect separator
            
            # Get sequence column (assume first column is sequences)
            seq_col = count_df.columns[0]
            sequences = count_df[seq_col].dropna().astype(str).tolist()
            
            # Take sample from this file
            all_sequences.extend(sequences[:50])  # 50 sequences per file
            files_sampled += 1
            
        except Exception as e:
            continue  # Skip problematic files
    
    if not all_sequences:
        raise ValueError("Could not read any sequences from count files")
    
    # Detect format
    detected_format = detect_sequence_format(all_sequences)
    
    if detected_format == 'mixed':
        raise ValueError("Mixed sequence formats detected in count files - ensure all files use consistent format")
    elif detected_format == 'unknown':
        raise ValueError("Could not determine sequence format from count files")
    elif detected_format not in ['dna', 'aa']:
        raise ValueError(f"Unexpected sequence format detected: {detected_format}")
    
    return detected_format
