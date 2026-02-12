"""
Unit tests for plotting functions in sortscore.visualization.plots.
"""
import pytest
import pandas as pd
import numpy as np
from sortscore.visualization.plots import plot_replicate_correlation

def test_plot_replicate_correlation(tmp_path):
    df = pd.DataFrame({
        'Rep1_score': [1, 2, 3, 4, 5],
        'Rep2_score': [1.1, 2.1, 2.9, 3.8, 5.2]
    })
    save_path = tmp_path / 'plot.png'
    # Should not raise
    plot_replicate_correlation(df, 'Rep1_score', 'Rep2_score', export=True, output=str(save_path))
    assert save_path.exists()

def test_plot_replicate_correlation_missing_column():
    df = pd.DataFrame({'other': [1, 2, 3]})
    with pytest.raises(ValueError):
        plot_replicate_correlation(df, x_col='Rep1_score', y_col='Rep2_score')
