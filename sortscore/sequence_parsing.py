"""
Sequence parsing utilities for sortscore package.

This module provides functions for DNA to protein translation and sequence comparison.

Examples
--------
>>> from sortscore.sequence_parsing import translate_dna, compare_to_reference
>>> translate_dna('ATGGCC')
'MA'
"""
from typing import List
from Bio.Seq import Seq

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

    Examples
    --------
    >>> compare_to_reference('MA', 'MT')
    'A.2.T'
    """
    differences = []
    min_length = min(len(ref_seq), len(sequence))
    for i in range(min_length):
        if ref_seq[i] != sequence[i]:
            difference = f'{ref_seq[i]}.{i+1}.{sequence[i]}'
            differences.append(difference)
    return ', '.join(differences)
