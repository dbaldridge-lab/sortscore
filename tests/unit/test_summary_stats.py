import pandas as pd

from sortscore.analysis.summary_stats import calculate_summary_stats


def test_calculate_summary_stats_returns_flat_sections():
    scores_df = pd.DataFrame(
        {
            "aa_seq_diff": ["=", "A.2.S", "Q.3.*"],
            "annotate_aa": ["synonymous", "missense_aa", "nonsense"],
            "annotate_dna": ["synonymous", "missense_dna", "missense_dna"],
            "score": [10.0, 20.0, 5.0],
        }
    )

    stats = calculate_summary_stats(
        scores_df,
        "score",
    )

    assert stats["overall"] == {"avg": 11.67, "median": 10.0, "min": 5.0, "max": 20.0, "std": 7.64}
    assert stats["synonymous_wt"] == {"avg": 10.0, "median": 10.0, "min": 10.0, "max": 10.0, "std": 0.0}
    assert stats["missense"] == {"avg": 20.0, "median": 20.0, "min": 20.0, "max": 20.0, "std": 0.0}
    assert stats["nonsense"] == {"avg": 5.0, "median": 5.0, "min": 5.0, "max": 5.0, "std": 0.0}
