from pathlib import Path

import pandas as pd
import pytest

from sortscore.analysis.aa_scores import build_aa_scores_table
from sortscore.analysis.batch_normalization import (
    _build_normalization_stats,
    run_batch_analysis,
    save_batch_results,
)


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


def test_build_normalization_stats_returns_nested_stage_and_final_sections():
    raw_scores = pd.DataFrame(
        {
            "batch": ["tile1", "tile1", "tile2", "tile2"],
            "aa_seq_diff": ["=", "Q.2.*", "=", "A.3.S"],
            "annotate_aa": ["synonymous", "nonsense", "synonymous", "missense_aa"],
            "annotate_dna": ["wt_dna", "missense_dna", "synonymous", "missense_dna"],
            "avgscore": [10.0, 2.0, 12.0, 20.0],
        }
    )
    final_scores = raw_scores.copy()
    final_scores["avgscore"] = [0.0, -2.0, 1.0, 3.0]

    stats = _build_normalization_stats(
        raw_tile_values={
            "tile1": {"wt_dna_score": 10.0, "non_avg": 2.0},
            "tile2": {"syn_median": 12.0},
        },
        raw_global_values={"wt_dna_score": 10.0},
        wt_stage_global_values={"syn_mean": 11.0, "syn_std": 1.4142135623730951},
        zscore_tile_values={
            "tile1": {"non_avg_zscore": -2.0},
            "tile2": {},
        },
        zscore_global_values={"non_avg_zscore": -2.0},
        final_scores=final_scores,
        wt_factors={"tile1": 1.1, "tile2": 0.9},
        path_factors={"tile1": 1.5, "tile2": 0.5},
    )

    assert "raw" in stats
    assert "wt_alignment" in stats
    assert "zscore" in stats
    assert stats["raw"]["global"]["wt_dna_score"] == 10.0
    assert stats["raw"]["tile1"]["wt_dna_score"] == 10.0
    assert stats["raw"]["tile2"]["syn_median"] == 12.0
    assert stats["wt_alignment"]["global"]["syn_mean"] == 11.0
    assert stats["wt_alignment"]["tile1"]["normalization_factor"] == 1.1
    assert stats["zscore"]["global"]["non_avg_zscore"] == -2.0
    assert stats["zscore"]["tile1"]["non_avg_zscore"] == -2.0
    assert stats["zscore"]["tile2"]["pathogenic_normalization_factor"] == 0.5
    assert stats["final"]["global"]["overall"] == {"avg": 0, "min": -2, "max": 3}
    assert stats["final"]["tile1"]["nonsense"] == {"avg": -2, "min": -2, "max": -2}


def test_run_batch_analysis_loads_batch_config_entries(tmp_path):
    tile_output_dir = tmp_path / "tile1"
    scores_dir = tile_output_dir / "scores"
    scores_dir.mkdir(parents=True)

    pd.DataFrame(
        {
            "variant_seq": ["AAA", "AAG", "TAA"],
            "dna_seq_diff": ["=", "A.1.G", "T.1.A"],
            "aa_seq_diff": ["=", "K.1.R", "Q.1.*"],
            "annotate_dna": ["wt_dna", "synonymous", "missense_dna"],
            "annotate_aa": ["synonymous", "synonymous", "nonsense"],
            "avgscore": [10.0, 12.0, 2.0],
            "avgscore_rep_weighted": [10.0, 12.0, 2.0],
            "Rep1.score": [10.0, 12.0, 2.0],
            "Rep2.score": [10.0, 12.0, 2.0],
        }
    ).to_csv(scores_dir / "tile1_dna_scores.csv", index=False)

    results = run_batch_analysis(
        {
            "batch_normalization_method": "zscore_2pole",
            "pathogenic_control_type": "nonsense",
            "combined_output_dir": str(tmp_path),
            "experiments": [
                {
                    "tile": 1,
                    "output_dir": str(tile_output_dir),
                    "avg_method": "simple-avg",
                    "min_pos": 1,
                    "max_pos": 3,
                    "wt_seq": "AAA",
                    "mutagenesis_type": "codon",
                }
            ],
        }
    )

    assert results["method"] == "zscore_2pole"
    assert results["output_dir"] == str((tmp_path / "normalized" / "zscore_2pole").resolve())
    assert len(results["experiments"]) == 1
    assert results["experiments"][0].min_pos == 1
    assert not results["normalized_scores"].empty
