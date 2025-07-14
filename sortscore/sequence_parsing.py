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
