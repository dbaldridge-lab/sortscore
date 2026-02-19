"""
Sequence parsing utilities for sortscore package.

This module provides functions for DNA to protein translation and sequence comparison.

Examples
--------
>>> from sortscore.sequence_parsing import translate_dna, compare_to_reference
>>> translate_dna('ATGGCC')
'MA'
"""
import re
from typing import List, Literal, Optional
from Bio.Seq import Seq


def _classify_single_sequence(seq: str) -> Optional[Literal['dna', 'aa']]:
    """Classify one sequence string as DNA, AA, or unknown."""
    seq = seq.strip()

    seq = seq.upper()

    aa_single_patterns = [
        r'^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY\*]$',
        r'^[ACDEFGHIKLMNPQRSTVWY]\.\d+\.[ACDEFGHIKLMNPQRSTVWY\*]$',
        r'^[ACDEFGHIKLMNPQRSTVWY]-\d+-[ACDEFGHIKLMNPQRSTVWY\*]$',
        r'^[ACDEFGHIKLMNPQRSTVWY]_\d+_[ACDEFGHIKLMNPQRSTVWY\*]$',
    ]
    for pattern in aa_single_patterns:
        if re.match(pattern, seq):
            return 'aa'

    aa_three_letter = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'TER'
    }
    match = re.match(r'^([A-Z]{3})(\d+)([A-Z]{3})$', seq)
    if match:
        from_aa, _, to_aa = match.groups()
        if from_aa in aa_three_letter and to_aa in aa_three_letter:
            return 'aa'

    dna_change_patterns = [
        r'^[ATCG]\d+[ATCG]$',
        r'^[ATCG]\.\d+\.[ATCG]$',
        r'^[ATCG]-\d+-[ATCG]$',
        r'^[ATCG]_\d+_[ATCG]$',
    ]
    for pattern in dna_change_patterns:
        if re.match(pattern, seq):
            return 'dna'

    seq_chars = set(seq)
    if seq_chars.issubset(set('ATCG')) and len(seq) > 1:
        return 'dna'
    if seq_chars.issubset(set('ACDEFGHIKLMNPQRSTVWY*')) and len(seq) > 1:
        return 'aa'
    if seq_chars.issubset(set('ATCGRYSWKMBDHVN')) and len(seq) > 1:
        return 'dna'
    return None


def detect_sequence_format(sequences: List[str]) -> Literal['dna', 'aa']:
    """
    Classify a set of sequences as DNA or amino-acid format.
    """
    if not sequences:
        raise ValueError("No sequences provided for format detection")

    valid_sequences = [seq for seq in sequences if seq and str(seq).strip()]
    if not valid_sequences:
        raise ValueError("No valid sequences provided for format detection")

    dna_count = 0
    aa_count = 0
    for seq in valid_sequences:
        fmt = _classify_single_sequence(str(seq))
        if fmt == 'dna':
            dna_count += 1
        elif fmt == 'aa':
            aa_count += 1

    total_classified = dna_count + aa_count
    if total_classified == 0:
        raise ValueError("All input sequences are unrecognized. Could not determine sequence format.")
    if dna_count == total_classified:
        return 'dna'
    if aa_count == total_classified:
        return 'aa'
    raise ValueError("Mixed sequence formats detected (both DNA and AA present). Ensure all sequences use a consistent format.")


def detect_sequence_format_from_counts(
    counts: dict,
    sample_size: int = 200,
    sequence_column: str = 'variant_seq',
) -> Literal['dna', 'aa']:
    """
    Detect sequence format from loaded count tables.

    Uses ``sequence_column`` when present, otherwise falls back to the first
    column in each table.
    """
    sampled_sequences = []
    for rep_dict in counts.values():
        for df in rep_dict.values():
            if df is None or df.empty:
                continue

            if sequence_column in df.columns:
                seq_col = sequence_column
            else:
                seq_col = df.columns[0]

            sampled_sequences.extend(df[seq_col].dropna().astype(str).head(sample_size).tolist())
            if len(sampled_sequences) >= sample_size:
                break
        if len(sampled_sequences) >= sample_size:
            break

    if not sampled_sequences:
        raise ValueError("No sequences available to detect sequence format from loaded counts.")

    return detect_sequence_format(sampled_sequences)


def translate_dna(dna_sequence: str) -> str:
    """
    Translate a DNA sequence into a protein sequence using standard genetic code.

    Parameters
    ----------
    dna_sequence : str
        DNA sequence (string of A, T, C, G).

    Returns
    -------
    protein_sequence : str
        Translated protein sequence.

    Examples
    --------
    >>> translate_dna('ATGGCC')
    'MA'
    """
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

def compare_to_reference(ref_seq: str, sequence: str) -> str:
    """
    Compare a protein sequence to a reference and report differences.
    
    Stop codons are represented as '*'.

    Parameters
    ----------
    ref_seq : str
        Reference protein sequence.
    sequence : str
        Protein sequence to compare.

    Returns
    -------
    differences : str
        Comma-separated string of differences in the format 'ref.position.alt'.
        Stop codons are shown as '*' (e.g., 'Q.10.*' for a stop-gained variant).

    Examples
    --------
    >>> compare_to_reference('MA', 'MT')
    'A.2.T'
    >>> compare_to_reference('MQ', 'MX')  # X (stop codon) becomes *
    'Q.2.*'
    """
    differences = []
    min_length = min(len(ref_seq), len(sequence))
    for i in range(min_length):
        if ref_seq[i] != sequence[i]:
            # Map X (stop codon) to * for standard notation
            ref_aa = '*' if ref_seq[i] == 'X' else ref_seq[i]
            var_aa = '*' if sequence[i] == 'X' else sequence[i]
            difference = f'{ref_aa}.{i+1}.{var_aa}'
            differences.append(difference)
    return ', '.join(differences)

def get_codons(dna_seq: str) -> List[str]:
    """
    Split a DNA sequence into codons (groups of 3 nucleotides).

    Parameters
    ----------
    dna_seq : str
        DNA sequence.

    Returns
    -------
    codons : List[str]
        List of codon strings.

    Raises
    ------
    ValueError
        If DNA sequence length is not divisible by 3.

    Examples
    --------
    >>> get_codons('ATGGCC')
    ['ATG', 'GCC']
    """
    if len(dna_seq) % 3 != 0:
        raise ValueError("The length of the DNA sequence must be divisible by 3 to label codons.")
    return [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]

def compare_codon_lists(wt_seq: str, variant_seq: str) -> str:
    """
    Compare codon sequences between wild-type and variant DNA sequences.

    Parameters
    ----------
    wt_seq : str
        Wild-type DNA sequence.
    variant_seq : str
        Variant DNA sequence.

    Returns
    -------
    differences : str
        Period-separated string of codon differences in format 'wtAA(wtCodon).position.varAA(varCodon)'.

    Examples
    --------
    >>> compare_codon_lists('ATGGCC', 'ATGTCC')
    'A(GCC).2.S(TCC)'
    """
    wt_codons = get_codons(wt_seq)
    variant_codons = get_codons(variant_seq)
    differences = []
    for i, (wt_codon, var_codon) in enumerate(zip(wt_codons, variant_codons)):
        if wt_codon != var_codon:
            wt_aa = translate_dna(wt_codon)
            var_aa = translate_dna(var_codon)
            differences.append(f"{wt_aa}({wt_codon}).{i+1}.{var_aa}({var_codon})")
    return '.'.join(differences)

def convert_aa_to_three_letter(aa_code: str) -> str:
    """
    Convert single-letter amino acid code to three-letter code.
    
    Parameters
    ----------
    aa_code : str
        Single-letter amino acid code (e.g., 'A', 'R', '*')
        
    Returns
    -------
    str
        Three-letter amino acid code (e.g., 'Ala', 'Arg', 'Ter')
        
    Examples
    --------
    >>> convert_aa_to_three_letter('A')
    'Ala'
    >>> convert_aa_to_three_letter('*')
    'Ter'
    """
    conversion_map = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
        '*': 'Ter'
    }
    return conversion_map.get(aa_code, aa_code)

def get_reference_sequence(wt_seq: str, target_format: str) -> str:
    """
    Get the appropriate reference sequence for comparison based on target format.
    
    Parameters
    ----------
    wt_seq : str
        Wild-type sequence (nucleotide or protein sequence)
    target_format : str
        Target format for analysis ('dna' or 'aa')
        
    Returns
    -------
    str
        Reference sequence in the target format
        
    Raises
    ------
    ValueError
        If conversion is invalid (e.g., trying to convert AA to DNA)
        
    Examples
    --------
    >>> get_reference_sequence('ATGGCC', 'aa')
    'MA'
    >>> get_reference_sequence('MKVL', 'aa')
    'MKVL'
    """
    wt_seq_format = detect_sequence_format([wt_seq])
    
    if target_format == 'aa':
        if wt_seq_format == 'dna':
            return translate_dna(wt_seq)
        else:
            return wt_seq
    elif target_format == 'dna':
        if wt_seq_format == 'dna':
            return wt_seq
        else:
            raise ValueError("DNA analysis requires DNA wt_seq, invalid sequence provided")
    else:
        raise ValueError(f"Unknown target format: {target_format}. Must be 'dna' or 'aa'.")
