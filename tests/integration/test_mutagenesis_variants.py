"""
Test mutagenesis_variants functionality in heatmap generation.
"""
import pandas as pd
import pytest
from sortscore.visualization.heatmap_matrix import dms_matrix_template, get_dropout


def test_dms_matrix_template_default():
    """Test that default template uses standard amino acids."""
    matrix = dms_matrix_template(10, 'aa')
    expected_rows = ['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*']
    assert list(matrix.index) == expected_rows
    assert len(matrix.columns) == 10


def test_dms_matrix_template_custom_variants():
    """Test that custom mutagenesis_variants affects matrix rows."""
    custom_variants = ['M', 'V', 'L', 'I', 'A']
    matrix = dms_matrix_template(10, 'aa', custom_variants)
    assert list(matrix.index) == custom_variants
    assert len(matrix.columns) == 10


def test_get_dropout_default():
    """Test dropout calculation with default variants."""
    # Create a 3x21 matrix (3 positions, 21 standard variants)
    data = pd.DataFrame(index=['W', 'F', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'Q', 'N', 'D', 'E', 'H', 'R', 'K', '*'],
                       columns=[1, 2, 3])
    data.iloc[0, 0] = 1.0  # Fill one cell
    
    dropout_num, dropout_percent = get_dropout(data)
    expected_total = 3 * 21  # 63 total possible variants
    expected_dropout = 62  # 62 missing variants
    
    assert dropout_num == expected_dropout
    assert dropout_percent == round((expected_dropout / expected_total) * 100, 1)


def test_get_dropout_custom_variants():
    """Test dropout calculation with custom variants."""
    custom_variants = ['M', 'V', 'L']
    data = pd.DataFrame(index=custom_variants, columns=[1, 2, 3])
    data.iloc[0, 0] = 1.0  # Fill one cell
    
    dropout_num, dropout_percent = get_dropout(data, custom_variants)
    expected_total = 3 * 3  # 9 total possible variants
    expected_dropout = 8  # 8 missing variants
    
    assert dropout_num == expected_dropout
    assert dropout_percent == round((expected_dropout / expected_total) * 100, 1)


def test_mutagenesis_variants_filtering():
    """Test that only specified variants appear in matrix."""
    hydrophobic_variants = ['M', 'I', 'L', 'V', 'F', 'W', 'Y', 'A']
    matrix = dms_matrix_template(5, 'aa', hydrophobic_variants)
    
    # Matrix should only have hydrophobic variants
    assert len(matrix.index) == len(hydrophobic_variants)
    assert all(aa in hydrophobic_variants for aa in matrix.index)
    
    # Should not have polar/charged variants
    assert 'K' not in matrix.index
    assert 'R' not in matrix.index
    assert 'D' not in matrix.index
    assert 'E' not in matrix.index