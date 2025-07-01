"""
Unit tests for all plotting functions in sortscore.visualization.plots.
"""
import pytest
import pandas as pd
import numpy as np
from sortscore.visualization import plots

def test_plot_activity_score_distribution(tmp_path):
    df = pd.DataFrame({'avgscore': np.random.normal(1, 0.2, 100)})
    save_path = tmp_path / 'dist.png'
    plots.plot_activity_score_distribution(df, score_col='avgscore', save_path=str(save_path))
    assert save_path.exists()

def test_plot_beeswarm(tmp_path):
    df = pd.DataFrame({
        'annotate_aa': ['A', 'B', 'A', 'B'],
        'score': [1.0, 2.0, 1.5, 2.5]
    })
    save_path = tmp_path / 'beeswarm.png'
    plots.plot_beeswarm(df, x='annotate_aa', y='score', export=True, output=str(save_path))
    assert save_path.exists()

def test_plot_histogram(tmp_path):
    df = pd.DataFrame({
        'annotate_dna': ['wt_dna', 'spike-in', 'other', 'other'],
        'score': [1.0, 0.5, 0.8, 1.2]
    })
    save_path = tmp_path / 'hist.png'
    plots.plot_histogram(df, score_col='score', group_col='annotate_dna', export=True, output=str(save_path))
    assert save_path.exists()

def test_plot_heatmap(tmp_path):
    # Minimal test for heatmap plotting
    df = pd.DataFrame({
        'aa_seq_diff': ['A.1.B', 'B.2.A'],
        'score': [1.0, 2.0]
    })
    save_path = tmp_path / 'heatmap.png'
    plots.plot_heatmap(
        data=df,
        score_col='score',
        num_aa=2,
        wt_seq='AB',
        min_pos=1,
        mutant_type='aa',
        wt_score=1.0,
        fig_size='small',
        export=True,
        output=str(save_path),
        tick_values=[0.5, 1.0, 2.0],
        tick_labels=['Low', 'WT', 'High'],
        motif_indices=None,
        row_avg=False,
        title='Test Heatmap'
    )
    assert save_path.exists()
