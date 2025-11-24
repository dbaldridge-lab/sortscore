"""
Sequence annotation utilities for Sort-seq variant analysis.

This module provides functions for annotating variant DataFrames with sequence differences,
translations, and other derived sequence information.

Examples
--------
>>> from sortscore.analysis.annotation import annotate_scores_dataframe
>>> annotated_df = annotate_scores_dataframe(scores_df, experiment)
"""
import pandas as pd
from typing import Optional
from sortscore.sequence_parsing import compare_to_reference, compare_codon_lists, translate_dna
import re


def add_hgvs_notations(scores_df: pd.DataFrame, wt_dna_seq: str, variant_type: str, transcript_id: Optional[str] = None) -> pd.DataFrame:
    """
    Add HGVS notation columns to scores DataFrame.
    
    For AA variant type: adds p. notation only
    For DNA variant type: adds both p. and c. notation
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with variant sequences and aa_seq_diff column
    wt_dna_seq : str
        Wild-type DNA sequence
    variant_type : str
        Either 'aa' or 'dna'
    transcript_id : str, optional
        Transcript ID for HGVS notation
        
    Returns
    -------
    pd.DataFrame
        DataFrame with added HGVS notation columns
    """
    df = scores_df.copy()
    
    if 'aa_seq_diff' not in df.columns:
        return df
    
    # Always add p. notation for protein changes
    def format_p_notation(aa_change):
        if not aa_change or aa_change == '':
            return ''
        
        # Convert "P.23.L" format to "p.Pro23Leu" format
        match = re.match(r'([A-Z*])\.(\d+)\.([A-Z*])', aa_change)
        if not match:
            return ''
        
        ref_aa, pos, alt_aa = match.groups()
        
        # Convert single letter to three letter code
        aa_map = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
            'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
            'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
            '*': 'Ter'
        }
        
        ref_long = aa_map.get(ref_aa, ref_aa)
        alt_long = aa_map.get(alt_aa, alt_aa)
        
        if ref_aa == alt_aa:
            p_notation = f"p.{ref_long}{pos}="
        else:
            p_notation = f"p.{ref_long}{pos}{alt_long}"
        
        if transcript_id:
            return f"{transcript_id}:{p_notation}"
        else:
            return p_notation
    
    df['hgvs_p'] = df['aa_seq_diff'].apply(format_p_notation)
    return df


def annotate_scores_dataframe(
    scores_df: pd.DataFrame, 
    wt_dna_seq: str, 
    variant_type: str = 'dna',
    transcript_id: Optional[str] = None
) -> pd.DataFrame:
    """
    Add sequence annotation columns to a scores DataFrame.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with variant sequences and scores.
    wt_dna_seq : str
        Wild-type DNA reference sequence.
    variant_type : str, default 'dna'
        Type of variants ('dna' or 'aa').
    transcript_id : str, optional
        Transcript ID for HGVS notation.
    
    Returns
    -------
    annotated_df : pd.DataFrame
        DataFrame with added annotation columns.
        
    Examples
    --------
    >>> annotated_df = annotate_scores_dataframe(scores_df, wt_seq, 'dna')
    """
    df = scores_df.copy()
    
    # Check if aa_seq_diff already exists (from pre-annotated data)
    has_pre_annotated_aa = 'aa_seq_diff' in df.columns
    
    if variant_type == 'dna':
        # Add codon differences
        df['codon_diff'] = df['variant_seq'].apply(
            lambda x: compare_codon_lists(wt_dna_seq, x)
        )
        df['codon_diff'] = df['codon_diff'].fillna('')
        
        # Add DNA sequence differences
        df['dna_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_dna_seq, x)
        )
        df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
        
        # Add AA sequence annotations only if not pre-annotated
        if not has_pre_annotated_aa:
            wt_aa_seq = translate_dna(wt_dna_seq)
            df['aa_seq'] = df['variant_seq'].apply(translate_dna)
            df['aa_seq_diff'] = df['aa_seq'].apply(
                lambda x: compare_to_reference(wt_aa_seq, x)
            )
            df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        
    elif variant_type == 'aa':
        # For AA variants, add sequence differences only if not pre-annotated
        if not has_pre_annotated_aa:
            wt_aa_seq = translate_dna(wt_dna_seq) if len(wt_dna_seq) % 3 == 0 else wt_dna_seq
            df['aa_seq_diff'] = df['variant_seq'].apply(
                lambda x: compare_to_reference(wt_aa_seq, x)
            )
            df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
    
    # Map stop codon representations to * for standard notation in aa_seq_diff column
    if 'aa_seq_diff' in df.columns:
        df['aa_seq_diff'] = df['aa_seq_diff'].str.replace('X', '*', regex=False)
        df['aa_seq_diff'] = df['aa_seq_diff'].str.replace('Ter', '*', regex=False)
    
    # Add functional annotations
    df = add_variant_categories(df)
    
    # Add HGVS notations
    df = add_hgvs_notations(df, wt_dna_seq, variant_type, transcript_id)
    
    return df


def add_sequence_differences(df: pd.DataFrame, wt_dna_seq: str, variant_type: str = 'dna') -> pd.DataFrame:
    """
    Add sequence difference columns to a DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with variant sequences.
    wt_dna_seq : str
        Wild-type DNA sequence.
    variant_type : str, default 'dna'
        Type of variants ('dna' or 'aa').
        
    Returns
    -------
    df : pd.DataFrame
        DataFrame with sequence difference columns added.
    """
    df = df.copy()
    
    if variant_type == 'dna':
        # Add DNA sequence differences
        df['dna_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_dna_seq, x)
        )
        df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
        
        # Add AA sequence differences
        wt_aa_seq = translate_dna(wt_dna_seq)
        df['aa_seq'] = df['variant_seq'].apply(translate_dna)
        df['aa_seq_diff'] = df['aa_seq'].apply(
            lambda x: compare_to_reference(wt_aa_seq, x)
        )
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        
    elif variant_type == 'aa':
        # For AA variants, sequences are already amino acids
        wt_aa_seq = translate_dna(wt_dna_seq) if len(wt_dna_seq) % 3 == 0 else wt_dna_seq
        df['aa_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_aa_seq, x)
        )
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
    
    return df


def add_variant_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add variant category annotations based on existing sequence difference columns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sequence difference columns (aa_seq_diff, dna_seq_diff).
        
    Returns
    -------
    df : pd.DataFrame
        DataFrame with variant category columns added.
    """
    df = df.copy()
    
    # Classify variants based on AA changes
    if 'aa_seq_diff' in df.columns:
        def classify_aa_variant(aa_diff, dna_diff=None):
            if not aa_diff or aa_diff == '':
                # Check if this is true WT (no DNA changes) or synonymous (DNA changes but same AA)
                if dna_diff is not None and (not dna_diff or dna_diff == ''):
                    return 'wt_dna'
                else:
                    return 'synonymous'
            elif '*' in aa_diff:
                return 'nonsense'
            else:
                return 'missense_aa'
        
        if 'dna_seq_diff' in df.columns:
            df['annotate_aa'] = df.apply(lambda row: classify_aa_variant(row['aa_seq_diff'], row['dna_seq_diff']), axis=1)
        else:
            df['annotate_aa'] = df['aa_seq_diff'].apply(classify_aa_variant)
    
    # Classify DNA variants 
    if 'dna_seq_diff' in df.columns:
        def classify_dna_variant(dna_diff, aa_diff):
            if not dna_diff or dna_diff == '':
                return 'wt_dna'
            elif not aa_diff or aa_diff == '':
                return 'synonymous'
            else:
                return 'missense_dna'
        
        if 'aa_seq_diff' in df.columns:
            df['annotate_dna'] = df.apply(lambda row: classify_dna_variant(row['dna_seq_diff'], row['aa_seq_diff']), axis=1)
        else:
            df['annotate_dna'] = df['dna_seq_diff'].apply(lambda x: 'missense_dna' if x else 'wt_dna')
    
    return df