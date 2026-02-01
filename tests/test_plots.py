"""
Unit tests for plotting functions in sortscore.visualization.plots.
"""
import pytest
import pandas as pd
from sortscore.visualization.plots import plot_activity_score_distribution

def test_plot_activity_score_distribution(tmp_path):
    df = pd.DataFrame({'avgscore': [1, 2, 2, 3, 3, 3, 4, 4, 5]})
    save_path = tmp_path / 'plot.png'
    # Should not raise
    plot_activity_score_distribution(df, score_col='avgscore', save_path=str(save_path))
    assert save_path.exists()

def test_plot_activity_score_distribution_missing_column():
    df = pd.DataFrame({'other': [1, 2, 3]})
    with pytest.raises(ValueError):
        plot_activity_score_distribution(df, score_col='avgscore')
