"""
Matrix and heatmap utilities for MAVE visualization.

This module provides functions to create and fill MAVE matrices for heatmap plotting, using sequence parsing and stats utilities.
"""
import pandas as pd
import numpy as np
from typing import List, Optional
from sortscore.sequence_parsing import translate_dna

def extract_position(sequence_diff: str):
    parts = sequence_diff.split('.')
    if len(parts) == 3 and parts[1].isdigit():
        return parts[2], int(parts[1])
    return None, None

def extract_value(cell: str):
    if cell == 'WT':
        return 'WT'
    elif pd.notna(cell):
        return cell
    else:
        return np.nan

def dms_matrix_template(num_aa: int, variant_type: str = 'aa') -> pd.DataFrame:
    # Always use 1-based indexing internally for matrix columns
    column_values = list(range(1, num_aa + 1))
    if variant_type == 'aa':
        row_labels = ['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']
    elif variant_type == 'dna':
        row_labels = [
            'W(TGG)', 'F(TTT)', 'F(TTC)', 'Y(TAT)', 'Y(TAC)', 'P(CCT)', 'P(CCC)', 'P(CCA)', 'P(CCG)', 'M(ATG)',
            'I(ATT)', 'I(ATC)', 'I(ATA)', 'L(TTA)', 'L(TTG)', 'L(CTT)', 'L(CTC)', 'L(CTA)', 'L(CTG)',
            'V(GTT)', 'V(GTC)', 'V(GTA)', 'V(GTG)', 'A(GCT)', 'A(GCC)', 'A(GCA)', 'A(GCG)', 'G(GGT)', 'G(GGC)',
            'G(GGA)', 'G(GGG)', 'C(TGT)', 'C(TGC)', 'S(TCT)', 'S(TCC)', 'S(TCA)', 'S(TCG)', 'S(AGT)', 'S(AGC)',
            'T(ACT)', 'T(ACC)', 'T(ACA)', 'T(ACG)', 'Q(CAA)', 'Q(CAG)', 'N(AAT)', 'N(AAC)', 'D(GAT)', 'D(GAC)',
            'E(GAA)', 'E(GAG)', 'H(CAT)', 'H(CAC)', 'R(CGT)', 'R(CGC)', 'R(CGA)', 'R(CGG)', 'R(AGA)', 'R(AGG)',
            'K(AAA)', 'K(AAG)', '*(TAA)', '*(TAG)', '*(TGA)'
        ]
    return pd.DataFrame(index=row_labels, columns=column_values)

def make_dms_matrix(
    data: pd.DataFrame,
    score_col: str,
    num_aa: int,
    wt_seq: str,
    variant_type: str = 'aa'
) -> pd.DataFrame:
    data = data.dropna(subset=[score_col])
    matrix = dms_matrix_template(num_aa, variant_type)
    diff_col = 'aa_seq_diff' if variant_type == 'aa' else 'codon_diff'
    for _, row in data.iterrows():
        char, col = extract_position(row[diff_col])
        if char in matrix.index and col in matrix.columns:
            matrix.at[char, col] = row[score_col]
    if variant_type == 'aa':
        for index, amino_acid in enumerate(wt_seq, start=1):
            if amino_acid in matrix.index and index in matrix.columns:
                matrix.at[amino_acid, index] = 'WT'
    if variant_type == 'dna':
        codon_enumeration = [((i // 3)+1, wt_seq[i:i+3]) for i in range(0, len(wt_seq), 3)]
        for index, codon in codon_enumeration:
            aa = translate_dna(codon)
            codon_index = f'{aa}({codon})'
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

def get_dropout(df: pd.DataFrame):
    """
    Calculate the number and percent of missing (dropout) variants in a MAVE matrix.

    Parameters
    ----------
    df : pandas.DataFrame
        MAVE matrix DataFrame.

    Returns
    -------
    dropout_num : int
        Number of missing variants.
    dropout_percent : float
        Percent of missing variants.
    """
    num_aa = df.shape[1]
    possible_sequences = num_aa * 20 + 1
    dropout_num = df.isna().sum().sum()
    dropout_percent = round((dropout_num / possible_sequences) * 100, 1)
    return dropout_num, dropout_percent
