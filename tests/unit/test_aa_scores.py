import pandas as pd
import pytest

from sortscore.analysis.aa_scores import build_aa_scores_table
from sortscore.analysis.variant_aggregation import aggregate_aa_data, aggregate_synonymous_variants


def assert_synonymous_no_diff_group_retained(df: pd.DataFrame, score_col: str = "avgscore") -> None:
    """Assert that synonymous no-diff rows collapse into a single retained AA group."""
    assert set(df["aa_seq_diff"]) == {"A.2.=", "A.2.S"}
    syn_rows = df[
        (df["aa_seq_diff"] == "A.2.=") &
        (df["annotate_aa"] == "synonymous")
    ]
    assert len(syn_rows) == 1
    syn_row = syn_rows.iloc[0]
    assert syn_row[score_col] == 12.0
    return syn_row


def test_build_aa_scores_table_retains_synonymous_no_diff_rows():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.2.=", "A.2.=", "A.2.S"],
            "annotate_aa": ["synonymous", "synonymous", "missense_aa"],
            "avgscore": [10.0, 14.0, 30.0],
            "avgscore_rep_weighted": [11.0, 15.0, 31.0],
            "count.r1b1": [100.0, 150.0, 300.0],
            "count.r1b2": [50.0, 75.0, 125.0],
        }
    )

    aa_scores = build_aa_scores_table(scores_df, "avgscore")

    # The two synonymous input rows should collapse into one retained no-diff AA row.
    syn_row = assert_synonymous_no_diff_group_retained(aa_scores)
    assert syn_row["avgscore_rep_weighted"] == 13.0
    assert syn_row["count.r1b1"] == 250.0
    assert syn_row["count.r1b2"] == 125.0


def test_aggregate_aa_data_retains_synonymous_no_diff_rows():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.2.=", "A.2.=", "A.2.S"],
            "annotate_aa": ["synonymous", "synonymous", "missense_aa"],
            "avgscore": [10.0, 14.0, 30.0],
        }
    )

    aa_data = aggregate_aa_data(scores_df, "avgscore")

    # AA aggregation should keep the synonymous no-diff group instead of dropping it.
    assert_synonymous_no_diff_group_retained(aa_data)


def test_aggregate_synonymous_variants_uses_explicit_column_rules():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.13.*", "A.13.*", "A.13.*"],
            "annotate_aa": ["nonsense", "nonsense", "nonsense"],
            "avgscore": [1.0, 2.0, 3.0],
            "Rep1.score": [2.0, 4.0, 6.0],
            "count.r1b1": [449.0, 45077.0, 1635.0],
            "count.r1b2": [10.0, 20.0, 30.0],
            "norm.count.r1b1": [4.0, 5.0, 6.0],
            "Rep1.sum": [7.0, 8.0, 9.0],
            "prop.r1b1": [0.1, 0.2, 0.3],
            "SD_rep": [1.0, 1.5, 2.0],
        }
    )

    aa_scores = aggregate_synonymous_variants(scores_df)
    row = aa_scores.iloc[0]

    assert row["avgscore"] == 2.0
    assert row["Rep1.score"] == 4.0
    assert row["count.r1b1"] == 47161.0
    assert row["count.r1b2"] == 60.0
    assert row["norm.count.r1b1"] == 15.0
    assert row["Rep1.sum"] == 24.0
    assert "prop.r1b1" not in aa_scores.columns
    assert "SD_rep" not in aa_scores.columns


def test_aggregate_synonymous_variants_rejects_columns_without_rules():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.13.*"],
            "annotate_aa": ["nonsense"],
            "unexpected_metric": [1.0],
        }
    )

    with pytest.raises(ValueError, match="unexpected_metric"):
        aggregate_synonymous_variants(scores_df)


def test_build_aa_scores_table_preserves_count_columns_for_aa_only_scores():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["A.2.S"],
            "annotate_aa": ["missense_aa"],
            "avgscore": [30.0],
            "avgscore_rep_weighted": [31.0],
            "count.r1b1": [300.0],
            "count.r1b2": [125.0],
        }
    )

    aa_scores = build_aa_scores_table(scores_df, "avgscore")
    row = aa_scores.iloc[0]

    assert row["count.r1b1"] == 300.0
    assert row["count.r1b2"] == 125.0
