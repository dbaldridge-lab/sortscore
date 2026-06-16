"""
Sequence annotation utilities for Sort-seq variant analysis.

This module provides functions for annotating variant DataFrames with sequence differences,
translations, and other derived sequence information.

Examples
--------
>>> from sortscore.analysis.annotation import annotate_scores_dataframe
>>> annotated_df = annotate_scores_dataframe(scores_df, experiment)
"""
from typing import Optional

import pandas as pd
from sortscore.utils.sequence_parsing import (
    compare_aa_from_dna_reference,
    compare_codon_lists,
    compare_to_reference,
    parse_sequence_difference,
    split_sequence_differences,
    translate_dna,
)

NO_DIFF_MARKER = '='


# TODO: #37 redundant, see if we can remove
def annotate_scores_dataframe(
    scores_df: pd.DataFrame, 
    wt_dna_seq: str, 
    mutagenesis_type: str = 'aa',
) -> pd.DataFrame:
    """
    Add sequence annotation columns to a scores DataFrame.
    
    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with variant sequences and scores.
    wt_dna_seq : str
        Wild-type DNA reference sequence.
    mutagenesis_type : str, default 'aa'
        Type of mutagenesis ('codon', 'snv', 'aa').
    
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
    
    # Treat 'dna' as a DNA-sequence variant type (full-length DNA sequences)
    if mutagenesis_type in {'codon', 'snv'}:
        # Add codon differences
        df['codon_diff'] = df['variant_seq'].apply(
            lambda x: compare_codon_lists(wt_dna_seq, x)
        )
        df['codon_diff'] = df['codon_diff'].fillna('')
        
        # Add DNA sequence differences
        df['dna_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_dna_seq, x)
        )
        
        # Add AA sequence annotations only if not pre-annotated
        if not has_pre_annotated_aa:
            df['aa_seq'] = df['variant_seq'].apply(translate_dna)
            df['aa_seq_diff'] = df['variant_seq'].apply(
                lambda x: compare_aa_from_dna_reference(wt_dna_seq, x)
            )
        
    elif mutagenesis_type == 'aa':
        # For AA variants, add sequence differences only if not pre-annotated
        if not has_pre_annotated_aa:
            wt_aa_seq = translate_dna(wt_dna_seq) if len(wt_dna_seq) % 3 == 0 else wt_dna_seq
            df['aa_seq_diff'] = df['variant_seq'].apply(
                lambda x: compare_to_reference(wt_aa_seq, x)
            )
    
    # Add functional annotations
    df = add_variant_categories(df)
    
    return df

# TODO: #37 isn't this redundant with similar functions
def add_sequence_differences(
    df: pd.DataFrame,
    wt_dna_seq: str,
    mutagenesis_type: str = 'aa',
) -> pd.DataFrame:
    """
    Add sequence difference columns to a DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with variant sequences.
    wt_dna_seq : str
        Wild-type DNA sequence.
    mutagenesis_type : str, default 'aa'
        Type of mutagenesis ('codon', 'snv', 'aa').
        
    Returns
    -------
    df : pd.DataFrame
        DataFrame with sequence difference columns added.
    """
    df = df.copy()
    
    if mutagenesis_type in {'codon', 'snv'}:
        # Add DNA sequence differences
        df['dna_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_dna_seq, x)
        )
        
        # Add AA sequence differences
        wt_aa_seq = translate_dna(wt_dna_seq)
        df['aa_seq'] = df['variant_seq'].apply(translate_dna)
        df['aa_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_aa_from_dna_reference(wt_dna_seq, x)
        )
        
    elif mutagenesis_type == 'aa':
        # For AA variants, sequences are already amino acids
        wt_aa_seq = translate_dna(wt_dna_seq) if len(wt_dna_seq) % 3 == 0 else wt_dna_seq
        df['aa_seq_diff'] = df['variant_seq'].apply(
            lambda x: compare_to_reference(wt_aa_seq, x)
        )
    
    return df

def classify_aa_variant(aa_diff: str, dna_diff: Optional[str] = None) -> str:
    """
    Classify an amino-acid-level sequence difference string.

    Parameters
    ----------
    aa_diff : str
        Amino-acid difference string. Accepted values are:
        - `'='` for true wild-type amino-acid sequence
        - one or more comma-separated changes in `ref.position.alt` format
          such as `'A.2.S'` or `'A.2.=, K.3.='`
    dna_diff : str, optional
        Matching DNA-level difference string for the same variant. When
        provided, exact `'='` distinguishes true wild type from synonymous
        amino-acid no-change strings.
    Returns
    -------
    str
        One of `'wt_dna'`, `'synonymous'`, `'nonsense'`, `'missense_aa'`,
        or `'multiple_aa'`.
    """
    aa_changes = split_sequence_differences(aa_diff)
    if not aa_changes:
        raise ValueError(
            f"Could not split amino acid changes from aa_diff '{aa_diff}'. "
            "Expected one or more changes separated by commas."
        )
    try:
        aa_no_change = aa_changes == [NO_DIFF_MARKER] or (
            all(parse_sequence_difference(change)[2] == NO_DIFF_MARKER for change in aa_changes)
        )
    except ValueError as exc:
        raise ValueError(
            f"Could not parse amino acid changes from aa_diff '{aa_diff}'. "
            "Each change must use 'ref.position.alt' format."
        ) from exc
    if len(aa_changes) > 1:
        return 'multiple_aa'
    if aa_no_change:
        # Check if this is true WT (no DNA changes) or synonymous (DNA changes but same AA)
        if dna_diff == NO_DIFF_MARKER:
            return 'wt_dna'
        else:
            return 'synonymous'
    elif '*' in aa_diff:
        return 'nonsense'
    else:
        return 'missense_aa'

def _classify_single_aa_change_dna_variant(dna_changes: list[str]) -> str:
    if len(dna_changes) == 1:
        return 'snv'
    if len(dna_changes) == 2:
        return 'dinucleotide'
    if len(dna_changes) == 3:
        return 'trinucleotide'
    raise ValueError(
        f"Expected 1-3 DNA changes for a single amino acid change, found {len(dna_changes)}."
    )


def classify_dna_variant(dna_diff: str, aa_diff: Optional[str] = None) -> str:
    dna_changes = split_sequence_differences(dna_diff)
    if not dna_changes:
        raise ValueError(
            f"Could not split DNA changes from dna_diff '{dna_diff}'. "
            "Expected one or more changes separated by commas."
        )
    dna_is_no_change = dna_changes == [NO_DIFF_MARKER]
    if NO_DIFF_MARKER in dna_changes and not dna_is_no_change:
        raise ValueError(
            f"Could not parse DNA changes from dna_diff '{dna_diff}'. "
            f"'{NO_DIFF_MARKER}' cannot be combined with other DNA changes."
        )
    if not dna_is_no_change:
        try:
            for change in dna_changes:
                parse_sequence_difference(change)
        except ValueError as exc:
            raise ValueError(
                f"Could not parse DNA changes from dna_diff '{dna_diff}'. "
                "Each change must use 'ref.position.alt' format."
            ) from exc
    aa_changes = split_sequence_differences(aa_diff)
    if not aa_changes:
        raise ValueError(
            f"Could not split amino acid changes from aa_diff '{aa_diff}'. "
            "Expected one or more changes separated by commas."
        )
    try:
        aa_no_change = aa_changes == [NO_DIFF_MARKER] or (
            all(parse_sequence_difference(change)[2] == NO_DIFF_MARKER for change in aa_changes)
        )
    except ValueError as exc:
        raise ValueError(
            f"Could not parse amino acid changes from aa_diff '{aa_diff}'. "
            "Each change must use 'ref.position.alt' format."
        ) from exc
    if len(aa_changes) > 1:
        return 'multiple_aa'
    if dna_is_no_change:
        return 'wt_dna'
    elif aa_no_change:
        return 'synonymous'
    else:
        return _classify_single_aa_change_dna_variant(dna_changes)


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
    
    # Classify DNA variants 
    if 'dna_seq_diff' in df.columns:
        if 'aa_seq_diff' in df.columns:
            df['annotate_dna'] = df.apply(lambda row: classify_dna_variant(row['dna_seq_diff'], row['aa_seq_diff']), axis=1)
        else:
            raise ValueError("dna_seq_diff requires aa_seq_diff for DNA variant classification")
    
    
    # Classify variants based on AA changes
    if 'aa_seq_diff' in df.columns:
        if 'dna_seq_diff' in df.columns:
            df['annotate_aa'] = df.apply(lambda row: classify_aa_variant(row['aa_seq_diff'], row['dna_seq_diff']), axis=1)
        else:
            df['annotate_aa'] = df['aa_seq_diff'].apply(classify_aa_variant)

    return df
