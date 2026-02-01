"""
Unit tests for score calculation functions in sortscore.analysis.score.
"""
import pytest
import pandas as pd
from sortscore.analysis.score import calculate_activity_scores

def test_calculate_activity_scores_simple_avg():
    df1 = pd.DataFrame({'variant': ['A', 'B'], 'count1': [10, 20]})
    df2 = pd.DataFrame({'variant': ['A', 'B'], 'count2': [30, 40]})
    # Merge columns to match expected input
    df1 = df1.rename(columns={'count1': 'sample1'})
    df2 = df2.rename(columns={'count2': 'sample2'})
    result = calculate_activity_scores([df1, df2], method='simple-avg')
    assert 'avgscore' in result.columns
    assert result.shape[0] == 2
    assert abs(result.loc[result['variant'] == 'A', 'avgscore'].values[0] - 20) < 1e-6
    assert abs(result.loc[result['variant'] == 'B', 'avgscore'].values[0] - 30) < 1e-6

def test_calculate_activity_scores_invalid_method():
    df = pd.DataFrame({'variant': ['A'], 'sample1': [1]})
    with pytest.raises(ValueError):
        calculate_activity_scores([df], method='invalid')
