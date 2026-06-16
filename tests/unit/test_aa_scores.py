import pandas as pd

from sortscore.analysis.aa_scores import build_aa_scores_table
from sortscore.analysis.variant_aggregation import aggregate_aa_data


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
        }
    )

    aa_scores = build_aa_scores_table(scores_df, "avgscore")

    # The two synonymous input rows should collapse into one retained no-diff AA row.
    syn_row = assert_synonymous_no_diff_group_retained(aa_scores)
    assert syn_row["avgscore_rep_weighted"] == 13.0


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
