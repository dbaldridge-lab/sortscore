"""
Matrix and heatmap utilities for MAVE visualization.

This module provides functions to create and fill MAVE matrices for heatmap plotting, using sequence parsing and stats utilities.
"""
import pandas as pd
import numpy as np
import re
from typing import List, Optional
from sortscore.utils.sequence_parsing import translate_dna, convert_aa_to_three_letter

def extract_position(sequence_diff: str):
    """Extract variant and position from standardized sequence difference format 'ref.position.alt'."""
    if pd.isna(sequence_diff):
        return None, None
    parts = sequence_diff.split('.')
    if len(parts) == 3 and parts[1].isdigit():
        return parts[2], int(parts[1])  # alt, position
    return None, None

def extract_value(cell: str):
    if cell == 'WT':
        return 'WT'
    elif pd.notna(cell):
        return cell
    else:
        return np.nan

def generate_codon_labels(three_letter_aa: bool = False) -> List[str]:
    """
    Generate default codon labels for DNA variant type.
    
    Parameters
    ----------
    three_letter_aa : bool, default False
        If True, use three-letter amino acid codes (e.g., 'Met(ATG)').
        If False, use single-letter codes (e.g., 'M(ATG)').
    
    Returns
    -------
    List[str]
        List of codon labels grouped by amino acid properties.
    """
    # Default codons grouped by amino acid properties
    codon_data = [
        ('M', 'ATG'), ('C', 'TGT'), ('C', 'TGC'), ('W', 'TGG'), 
        ('F', 'TTT'), ('F', 'TTC'), ('Y', 'TAT'), ('Y', 'TAC'), 
        ('P', 'CCT'), ('P', 'CCC'), ('P', 'CCA'), ('P', 'CCG'),
        ('I', 'ATT'), ('I', 'ATC'), ('I', 'ATA'), 
        ('L', 'TTA'), ('L', 'TTG'), ('L', 'CTT'), ('L', 'CTC'), ('L', 'CTA'), ('L', 'CTG'),
        ('V', 'GTT'), ('V', 'GTC'), ('V', 'GTA'), ('V', 'GTG'), 
        ('A', 'GCT'), ('A', 'GCC'), ('A', 'GCA'), ('A', 'GCG'), 
        ('G', 'GGT'), ('G', 'GGC'), ('G', 'GGA'), ('G', 'GGG'), 
        ('S', 'TCT'), ('S', 'TCC'), ('S', 'TCA'), ('S', 'TCG'), ('S', 'AGT'), ('S', 'AGC'),
        ('T', 'ACT'), ('T', 'ACC'), ('T', 'ACA'), ('T', 'ACG'), 
        ('Q', 'CAA'), ('Q', 'CAG'), ('N', 'AAT'), ('N', 'AAC'), 
        ('D', 'GAT'), ('D', 'GAC'), ('E', 'GAA'), ('E', 'GAG'), 
        ('H', 'CAT'), ('H', 'CAC'), 
        ('R', 'CGT'), ('R', 'CGC'), ('R', 'CGA'), ('R', 'CGG'), ('R', 'AGA'), ('R', 'AGG'),
        ('K', 'AAA'), ('K', 'AAG'), ('*', 'TAA'), ('*', 'TAG'), ('*', 'TGA')
    ]
    
    if three_letter_aa:
        return [f'{convert_aa_to_three_letter(aa)}({codon})' for aa, codon in codon_data]
    else:
        return [f'{aa}({codon})' for aa, codon in codon_data]

def dms_matrix_template(num_positions: int, mutagenesis_type: str = 'aa', mutagenesis_variants: list = None, three_letter_aa: bool = False) -> pd.DataFrame:
    column_values = list(range(1, num_positions + 1))
    
    # Row labels depend on mutagenesis_type and mutagenesis_variants
    if mutagenesis_type == 'aa':
        if mutagenesis_variants is not None:
            if three_letter_aa:
                row_labels = [convert_aa_to_three_letter(aa) for aa in mutagenesis_variants]
            else:
                row_labels = mutagenesis_variants
        else:
            default_aa = ['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']
            if three_letter_aa:
                row_labels = [convert_aa_to_three_letter(aa) for aa in default_aa]
            else:
                row_labels = default_aa
    elif mutagenesis_type == 'codon':
        row_labels = generate_codon_labels(three_letter_aa)
    elif mutagenesis_type in {'snv'}:
        row_labels = mutagenesis_variants if mutagenesis_variants is not None else ['A', 'C', 'G', 'T']
    return pd.DataFrame(index=row_labels, columns=column_values)

def make_dms_matrix(
    data: pd.DataFrame,
    score_col: str,
    num_positions: int,
    wt_seq: str,
    mutagenesis_type: str = 'aa',
    mutagenesis_variants: list = None,
    three_letter_aa: bool = False
) -> pd.DataFrame:
    
    data = data.dropna(subset=[score_col])
    matrix = dms_matrix_template(num_positions, mutagenesis_type, mutagenesis_variants, three_letter_aa)
    
    if mutagenesis_type == 'snv':
        diff_col = 'dna_seq_diff'
    elif mutagenesis_type == 'codon':
        diff_col = 'codon_diff'
    else:
        diff_col = 'aa_seq_diff'
    
    for _, row in data.iterrows():
        char, col = extract_position(row[diff_col])
        if char in matrix.index and col in matrix.columns:
            matrix.at[char, col] = row[score_col]
    # Mark WT positions
    if mutagenesis_type == 'snv':
        # For DNA positions, mark WT bases at each nucleotide position
        for index, base in enumerate(wt_seq, start=1):
            if base in matrix.index and index in matrix.columns:
                matrix.at[base, index] = 'WT'
    elif mutagenesis_type == 'aa':
        # For AA positions, mark WT amino acids
        # Check if wt_seq is DNA (length divisible by 3 and contains only ATCG) or already AA
        is_dna = len(wt_seq) % 3 == 0 and all(c in 'ATCG' for c in wt_seq.upper())
        wt_aa_seq = translate_dna(wt_seq) if is_dna else wt_seq
        for index, amino_acid in enumerate(wt_aa_seq, start=1):
            # Convert to three-letter if needed for matrix lookup
            aa_key = convert_aa_to_three_letter(amino_acid) if three_letter_aa else amino_acid
            if aa_key in matrix.index and index in matrix.columns:
                matrix.at[aa_key, index] = 'WT'
    elif mutagenesis_type == 'codon':
        # For DNA variant type with AA positions (codon-level)
        codon_enumeration = [((i // 3)+1, wt_seq[i:i+3]) for i in range(0, len(wt_seq), 3)]
        for index, codon in codon_enumeration:
            aa = translate_dna(codon)
            # Convert to three-letter if needed for codon label
            aa_display = convert_aa_to_three_letter(aa) if three_letter_aa else aa
            codon_index = f'{aa_display}({codon})'
            if codon_index in matrix.index and index in matrix.columns:
                matrix.at[codon_index, index] = 'WT'
    return matrix

def fill_wt(dms_matrix: pd.DataFrame, wt_dna_score: float) -> pd.DataFrame:
    pd.set_option('future.no_silent_downcasting', True)
    heatmap_df = dms_matrix.map(extract_value)
    heatmap_df = heatmap_df.replace('WT', wt_dna_score).astype(float)
    return heatmap_df

def make_col_avg_df(heatmap_df: pd.DataFrame) -> pd.DataFrame:
    col_avg = heatmap_df.iloc[:-1].mean()
    col_avg_df = pd.DataFrame(col_avg).transpose()
    return col_avg_df

def get_dropout(df: pd.DataFrame, mutagenesis_variants: list = None):
    """
    Calculate the number and percent of missing (dropout) variants in a MAVE matrix.

    Parameters
    ----------
    df : pandas.DataFrame
        MAVE matrix DataFrame.
    mutagenesis_variants : list, optional
        List of variant amino acids being studied. If None, uses default 20 AAs + stop.

    Returns
    -------
    dropout_num : int
        Number of missing variants.
    dropout_percent : float
        Percent of missing variants.
    """
    num_aa = df.shape[1]
    if mutagenesis_variants is not None:
        num_variants = len(mutagenesis_variants)
    else:
        num_variants = 21  # Default: 20 AAs + stop codon
    
    possible_sequences = num_aa * num_variants
    dropout_num = df.isna().sum().sum()
    dropout_percent = round((dropout_num / possible_sequences) * 100, 1)
    return dropout_num, dropout_percent
