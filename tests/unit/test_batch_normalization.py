from pathlib import Path

import pandas as pd
import pytest

from sortscore.analysis.aa_scores import build_aa_scores_table
from sortscore.analysis.batch_normalization import save_batch_results


def test_save_batch_results_preserves_decimal_score_precision(tmp_path):
    output_dir = Path(tmp_path) / "normalized" / "zscore_2pole"
    normalized_scores = pd.DataFrame(
        {
            "variant_seq": ["AAA"],
            "batch": ["tile1"],
            "score.r1b1": [1.25],
            "score.r1b2": [2.75],
            "Rep1.score": [4.125],
            "avgscore": [3.875],
            "avgscore_rep_weighted": [3.625],
        }
    )
    results = {
        "normalized_scores": normalized_scores,
        "normalized_aa_scores": pd.DataFrame(),
        "combined_stats": {"example": 1},
    }

    save_batch_results(results, str(output_dir))

    saved_scores = pd.read_csv(output_dir / "scores" / "batch_scores.csv")
    for col, expected in {
        "score.r1b1": 1.25,
        "score.r1b2": 2.75,
        "Rep1.score": 4.125,
        "avgscore": 3.875,
        "avgscore_rep_weighted": 3.625,
    }.items():
        assert saved_scores.loc[0, col] == expected


def test_build_aa_scores_table_preserves_decimal_stats_for_codon_aggregation():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.1.B", "A.1.B"],
            "annotate_aa": ["missense_aa", "missense_aa"],
            "avgscore": [10.0, 13.0],
            "avgscore_rep_weighted": [10.5, 13.5],
            "Rep1.score": [9.0, 12.0],
            "Rep2.score": [12.0, 15.0],
        }
    )

    aa_scores = build_aa_scores_table(scores_df, "avgscore")
    row = aa_scores.iloc[0]

    assert row["SD_codon"] == pytest.approx(2.1213203435596424)
    assert row["SD_rep"] == pytest.approx(2.1213203435596424)
    assert row["SEM"] == pytest.approx(1.5)
    assert row["CI_lower"] == pytest.approx(7.226330542073606)
    assert row["CI_upper"] == pytest.approx(16.773669457926394)


def test_build_aa_scores_table_preserves_decimal_stats_for_aa_only_scores():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.1.B"],
            "annotate_aa": ["missense_aa"],
            "avgscore": [10.5],
            "avgscore_rep_weighted": [10.25],
            "Rep1.score": [9.0],
            "Rep2.score": [12.0],
        }
    )

    aa_scores = build_aa_scores_table(scores_df, "avgscore")
    row = aa_scores.iloc[0]

    assert row["SD_rep"] == pytest.approx(2.1213203435596424)
    assert row["SEM"] == pytest.approx(1.5)
    assert row["CI_lower"] == pytest.approx(-8.559307104263047)
    assert row["CI_upper"] == pytest.approx(29.559307104263045)
