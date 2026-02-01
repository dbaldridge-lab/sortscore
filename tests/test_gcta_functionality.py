"""
Test GCTA (DNA-level) functionality in heatmap generation.
"""
import pandas as pd
import pytest
from sortscore.visualization.heatmap_matrix import dms_matrix_template, make_dms_matrix, get_dropout


def test_gcta_matrix_template():
    """Test that GCTA matrix template works for DNA positions."""
    gcta_variants = ['G', 'C', 'T', 'A']
    matrix = dms_matrix_template(100, 'dna', gcta_variants, 'dna')
    
    # Should have 4 rows (GCTA) and 100 columns (DNA positions)
    assert matrix.shape == (4, 100)
    assert list(matrix.index) == gcta_variants
    assert list(matrix.columns) == list(range(1, 101))


def test_dna_vs_aa_positions():
    """Test difference between DNA and AA position types."""
    variants = ['M', 'V', 'L']
    
    # AA positions
    aa_matrix = dms_matrix_template(10, 'aa', variants, 'aa')
    assert aa_matrix.shape == (3, 10)
    
    # DNA positions (same variants, different positions)
    dna_matrix = dms_matrix_template(30, 'aa', variants, 'dna')
    assert dna_matrix.shape == (3, 30)


def test_make_dms_matrix_dna_positions():
    """Test make_dms_matrix with DNA positions and differences."""
    # Create test data with DNA differences
    data = pd.DataFrame({
        'variant_seq': ['ATCG', 'TTCG', 'ATGG'],
        'dna_seq_diff': ['A.1.T', 'C.3.G', ''],  # DNA differences
        'score': [1.5, 2.0, 1.0]
    })
    
    gcta_variants = ['G', 'C', 'T', 'A']
    wt_seq = 'ATCG'
    
    matrix = make_dms_matrix(
        data, 'score', 4, wt_seq, 'dna', gcta_variants, 'dna'
    )
    
    # Check that scores are placed correctly
    assert matrix.at['T', 1] == 1.5  # A.1.T -> T at position 1
    assert matrix.at['G', 3] == 2.0  # C.3.G -> G at position 3
    
    # Check WT positions
    assert matrix.at['A', 1] == 'WT'  # WT A at position 1
    assert matrix.at['T', 2] == 'WT'  # WT T at position 2
    assert matrix.at['C', 3] == 'WT'  # WT C at position 3
    assert matrix.at['G', 4] == 'WT'  # WT G at position 4


def test_dropout_calculation_gcta():
    """Test dropout calculation for GCTA experiments."""
    gcta_variants = ['G', 'C', 'T', 'A']
    
    # Create a 4x10 matrix (4 bases, 10 positions)
    data = pd.DataFrame(index=gcta_variants, columns=list(range(1, 11)))
    data.iloc[0, 0] = 1.0  # Fill only one cell
    
    dropout_num, dropout_percent = get_dropout(data, gcta_variants)
    expected_total = 10 * 4  # 40 total possible variants
    expected_dropout = 39  # 39 missing variants
    
    assert dropout_num == expected_dropout
    assert dropout_percent == round((expected_dropout / expected_total) * 100, 1)


def test_custom_dna_variants():
    """Test custom DNA variants beyond GCTA."""
    custom_variants = ['G', 'A', 'T']  # Only 3 bases
    matrix = dms_matrix_template(5, 'dna', custom_variants, 'dna')
    
    assert matrix.shape == (3, 5)
    assert list(matrix.index) == custom_variants
    assert 'C' not in matrix.index


def test_mixed_position_variant_types():
    """Test different combinations of position_type and variant_type."""
    # DNA variants with DNA positions (GCTA)
    matrix1 = dms_matrix_template(10, 'dna', ['G', 'C', 'T', 'A'], 'dna')
    assert matrix1.shape == (4, 10)
    
    # AA variants with DNA positions (e.g., for deep scanning)
    matrix2 = dms_matrix_template(30, 'aa', ['M', 'V', 'L'], 'dna')
    assert matrix2.shape == (3, 30)
    
    # DNA variants with AA positions (traditional codon-level)
    matrix3 = dms_matrix_template(10, 'dna', None, 'aa')
    assert matrix3.shape[1] == 10  # Should have 10 AA positions
    assert matrix3.shape[0] > 20  # Should have many codon variants


def test_gcta_example_config():
    """Test realistic GCTA experiment configuration."""
    # Simulate GCTA experiment: single nucleotide changes across DNA sequence
    wt_dna = 'ATGCGTAAC'  # 9 nucleotides
    gcta_variants = ['G', 'C', 'T', 'A']
    
    # Test data with single nucleotide changes
    test_data = pd.DataFrame({
        'variant_seq': ['GTGCGTAAC', 'ATCCGTAAC', 'ATGCTTAAC'],
        'dna_seq_diff': ['A.1.G', 'G.4.C', 'G.5.T'],
        'score': [0.8, 1.2, 0.5]
    })
    
    matrix = make_dms_matrix(
        test_data, 'score', 9, wt_dna, 'dna', gcta_variants, 'dna'
    )
    
    # Verify matrix dimensions
    assert matrix.shape == (4, 9)
    
    # Verify scores placed correctly
    assert matrix.at['G', 1] == 0.8  # A.1.G
    assert matrix.at['C', 4] == 1.2  # G.4.C
    assert matrix.at['T', 5] == 0.5  # G.5.T
    
    # Verify WT positions
    assert matrix.at['A', 1] == 'WT'  # WT A at position 1
    assert matrix.at['T', 2] == 'WT'  # WT T at position 2
    assert matrix.at['G', 3] == 'WT'  # WT G at position 3