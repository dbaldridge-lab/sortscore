import pandas as pd

from sortscore.analysis.summary_stats import calculate_summary_stats


def test_calculate_summary_stats_returns_flat_sections():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["=", "A.2.S", "Q.3.*"],
            "annotate_aa": ["synonymous", "missense_aa", "nonsense"],
            "annotate_dna": ["synonymous", "missense_dna", "missense_dna"],
            "avgscore": [10.0, 20.0, 5.0],
            "avgscore_rep_weighted": [11.0, 21.0, 6.0],
        }
    )

    stats = calculate_summary_stats(
        scores_df,
        "avgscore",
    )

    assert stats["overall"] == {"avg": 12, "min": 5, "max": 20}
    assert stats["synonymous_wt"] == {"avg": 10, "min": 10, "max": 10}
    assert stats["missense"] == {"avg": 20, "min": 20, "max": 20}
    assert stats["nonsense"] == {"avg": 5, "min": 5, "max": 5}
