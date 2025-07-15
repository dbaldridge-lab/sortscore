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


def annotate_scores_dataframe(
    scores_df: pd.DataFrame, 
    wt_dna_seq: str, 
    variant_type: str = 'dna'
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
    
    Returns
    -------
    annotated_df : pd.DataFrame
        DataFrame with added annotation columns.
        
    Examples
    --------
    >>> annotated_df = annotate_scores_dataframe(scores_df, wt_seq, 'dna')
    """
    df = scores_df.copy()
    
    if variant_type == 'dna':
        # Add codon differences
        df['codon_diff'] = df['variant_seq'].apply(
            lambda x: compare_codon_lists(wt_dna_seq, x)
        )
        df['codon_diff'] = df['codon_diff'].fillna('')
        
        # Add AA sequence annotations
        wt_aa_seq = translate_dna(wt_dna_seq)
        df['aa_seq'] = df['variant_seq'].apply(translate_dna)
        df['aa_seq_diff'] = df['aa_seq'].apply(
            lambda x: compare_to_reference(wt_aa_seq, x)
        )
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        
        # Add DNA sequence differences
        df['dna_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_dna_seq, x)
        )
        df['dna_seq_diff'] = df['dna_seq_diff'].fillna('')
        
        # Add functional annotations
        df = basic_annotate(df, wt_dna_seq)
        
    elif variant_type == 'aa':
        # For AA variants, sequences are already amino acids
        wt_aa_seq = translate_dna(wt_dna_seq) if len(wt_dna_seq) % 3 == 0 else wt_dna_seq
        df['aa_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_aa_seq, x)
        )
        df['aa_seq_diff'] = df['aa_seq_diff'].fillna('')
        
        # Add functional annotations for AA variants
        df = basic_annotate(df, wt_dna_seq)
    
    return df


def basic_annotate(df: pd.DataFrame, wt_dna_seq: str) -> pd.DataFrame:
    """
    Add basic functional annotation categories.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sequence difference columns.
    wt_dna_seq : str
        Wild-type DNA sequence.
        
    Returns
    -------
    annotated_df : pd.DataFrame
        DataFrame with functional annotation columns.
    """
    df = df.copy()
    
    # Basic annotation logic
    if 'aa_seq_diff' in df.columns:
        # Classify variants based on AA changes
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
    
    if 'dna_seq_diff' in df.columns:
        # Classify DNA variants 
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


def aggregate_synonymous_variants(scores_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate synonymous variants by averaging their scores.
    
    Groups variants by their AA sequence difference and annotation type,
    then averages activity scores across synonymous variants to create
    a single representative row for each unique AA change.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with annotated scores including 'aa_seq_diff' and 'annotate_aa' columns.
        
    Returns
    -------
    aggregated_df : pd.DataFrame
        DataFrame with synonymous variants averaged into single rows.
        
    Examples
    --------
    >>> aggregated = aggregate_synonymous_variants(annotated_scores)
    """
    if 'aa_seq_diff' not in scores_df.columns or 'annotate_aa' not in scores_df.columns:
        raise ValueError("DataFrame must contain 'aa_seq_diff' and 'annotate_aa' columns")
    
    # Define columns to aggregate (all numeric columns)
    numeric_cols = scores_df.select_dtypes(include=['number']).columns.tolist()
    
    # Group by AA sequence difference and annotation type, then average
    aggregated_df = scores_df.groupby(['aa_seq_diff', 'annotate_aa'])[numeric_cols].mean().reset_index()
    
    # Add back non-numeric columns from the first occurrence (they should be the same for synonymous variants)
    other_cols = [col for col in scores_df.columns if col not in numeric_cols and 
                  col not in ['aa_seq_diff', 'annotate_aa']]
    if other_cols:
        first_others = scores_df.groupby(['aa_seq_diff', 'annotate_aa'])[other_cols].first().reset_index()
        aggregated_df = aggregated_df.merge(first_others, on=['aa_seq_diff', 'annotate_aa'])
    
    return aggregated_df


def add_functional_annotations(
    df: pd.DataFrame, 
    wt_dna_seq: str, 
    spike_in_seq: Optional[str] = None
) -> pd.DataFrame:
    """
    Add functional annotation categories (synonymous, missense, nonsense, etc.).
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sequence difference columns.
    wt_dna_seq : str
        Wild-type DNA sequence.
    spike_in_seq : str, optional
        Spike-in sequence for identification.
        
    Returns
    -------
    annotated_df : pd.DataFrame
        DataFrame with functional annotation columns.
    """
    return basic_annotate(df, wt_dna_seq)