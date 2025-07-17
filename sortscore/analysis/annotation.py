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
    
    # Add functional annotations
    df = add_variant_categories(df)
    
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


def aggregate_synonymous_variants(scores_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate synonymous variants by averaging their scores.
    
    Groups variants by their AA sequence difference and annotation type,
    then averages individual replicate scores and recalculates the three
    avgscore columns.
    
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
    import numpy as np
    
    if 'aa_seq_diff' not in scores_df.columns or 'annotate_aa' not in scores_df.columns:
        raise ValueError("DataFrame must contain 'aa_seq_diff' and 'annotate_aa' columns")
    
    # First, average the individual replicate scores and other columns (excluding avgscore columns)
    cols_to_average = [col for col in scores_df.columns if col not in ['aa_seq_diff', 'annotate_aa', 'avgscore', 'avgscore_rep_weighted', 'avgscore_codon_weighted']]
    numeric_cols = scores_df[cols_to_average].select_dtypes(include=['number']).columns.tolist()
    
    # Group by AA sequence difference and annotation type, then average
    aggregated_df = scores_df.groupby(['aa_seq_diff', 'annotate_aa'])[numeric_cols].mean().reset_index()
    
    # Recalculate the three avgscore columns from the averaged replicate scores
    # Simple average
    rep_score_cols = [col for col in aggregated_df.columns if 'Rep' in col and 'score' in col and 'cw' not in col and 'sum' not in col]
    if rep_score_cols:
        aggregated_df['avgscore'] = aggregated_df[rep_score_cols].mean(axis=1)
    
    # Rep-weighted average (need to recalculate weights from the aggregated rep sums)
    rep_sum_cols = [col for col in aggregated_df.columns if 'Rep' in col and 'sum' in col and 'syn' not in col]
    if rep_sum_cols and rep_score_cols:
        total_weight = aggregated_df[rep_sum_cols].sum(axis=1)
        weighted_sum = sum(aggregated_df[score_col] * aggregated_df[weight_col] 
                          for score_col, weight_col in zip(rep_score_cols, rep_sum_cols))
        aggregated_df['avgscore_rep_weighted'] = weighted_sum / total_weight
    
    # Codon-weighted average (use individual codon-weighted scores)
    rep_score_cw_cols = [col for col in aggregated_df.columns if 'Rep' in col and 'score.cw' in col]
    if rep_score_cw_cols:
        aggregated_df['avgscore_codon_weighted'] = aggregated_df[rep_score_cw_cols].sum(axis=1)
    
    # Add back non-numeric columns from the first occurrence
    other_cols = [col for col in scores_df.columns if col not in numeric_cols and 
                  col not in ['aa_seq_diff', 'annotate_aa']]
    if other_cols:
        first_others = scores_df.groupby(['aa_seq_diff', 'annotate_aa'])[other_cols].first().reset_index()
        aggregated_df = aggregated_df.merge(first_others, on=['aa_seq_diff', 'annotate_aa'])
    
    return aggregated_df