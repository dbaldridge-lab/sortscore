import numpy as np
import pandas as pd

from sortscore.analysis.score import calculate_full_activity_scores


def _scoring_inputs(replicates=(1,), bins=(1, 2, 3, 4)):
    counts = {}
    mfi = {}
    total_reads = {}

    for rep in replicates:
        counts[rep] = {}
        mfi[rep] = {}
        total_reads[rep] = {}

        for bin_number in bins:
            counts[rep][bin_number] = pd.DataFrame({"count": [10.0]})
            mfi[rep][bin_number] = float(bin_number)
            total_reads[rep][bin_number] = 1000.0

    return counts, mfi, total_reads


def test_zero_counts_do_not_satisfy_bins_required_or_assign_score():
    counts, mfi, total_reads = _scoring_inputs()
    merged = pd.DataFrame(
        {
            "variant_seq": ["variant2"],
            "count.r1b1": [10.0],
            "count.r1b2": [0.0],
            "count.r1b3": [0.0],
            "count.r1b4": [0.0],
        }
    )

    result = calculate_full_activity_scores(
        counts=counts,
        mfi=mfi,
        total_reads=total_reads,
        merged_df=merged,
        min_bins=3,
        min_reps=1,
        minread_threshold=0,
    )

    row = result.iloc[0]
    assert row["variant_seq"] == "variant2"
    assert row["count.r1b2"] == 0.0
    assert np.isnan(row["norm.count.r1b2"])
    assert np.isnan(row["Rep1.score"])
    assert np.isnan(row["score"])


def test_zero_count_bins_are_preserved_as_raw_counts_but_not_normalized_counts():
    counts, mfi, total_reads = _scoring_inputs()
    merged = pd.DataFrame(
        {
            "variant_seq": ["variant1"],
            "count.r1b1": [10.0],
            "count.r1b2": [20.0],
            "count.r1b3": [0.0],
            "count.r1b4": [30.0],
        }
    )

    result = calculate_full_activity_scores(
        counts=counts,
        mfi=mfi,
        total_reads=total_reads,
        merged_df=merged,
        min_bins=3,
        min_reps=1,
        minread_threshold=0,
    )

    row = result.iloc[0]
    assert row["count.r1b3"] == 0.0
    assert np.isnan(row["norm.count.r1b3"])
    assert pd.notna(row["Rep1.score"])


def test_minread_threshold_masks_normalized_values_without_changing_raw_counts():
    counts, mfi, total_reads = _scoring_inputs()
    total_reads = {1: {bin_number: 1000000 for bin_number in (1, 2, 3, 4)}}
    merged = pd.DataFrame(
        {
            "variant_seq": ["variant1"],
            "count.r1b1": [10.0],
            "count.r1b2": [1.0],
            "count.r1b3": [10.0],
            "count.r1b4": [10.0],
        }
    )

    result = calculate_full_activity_scores(
        counts=counts,
        mfi=mfi,
        total_reads=total_reads,
        merged_df=merged,
        min_bins=3,
        min_reps=1,
        minread_threshold=2,
    )

    row = result.iloc[0]
    assert row["count.r1b2"] == 1.0
    assert np.isnan(row["norm.count.r1b2"])
    assert pd.notna(row["Rep1.score"])


def test_replicates_that_fail_bins_required_do_not_satisfy_reps_required_or_assign_score():
    counts, mfi, total_reads = _scoring_inputs(replicates=(1, 2, 3), bins=(1, 2, 3))
    merged = pd.DataFrame(
        {
            "variant_seq": ["variant1"],
            "count.r1b1": [10.0],
            "count.r1b2": [0.0],
            "count.r1b3": [0.0],
            "count.r2b1": [0.0],
            "count.r2b2": [10.0],
            "count.r2b3": [0.0],
            "count.r3b1": [0.0],
            "count.r3b2": [0.0],
            "count.r3b3": [10.0],
        }
    )

    result = calculate_full_activity_scores(
        counts=counts,
        mfi=mfi,
        total_reads=total_reads,
        merged_df=merged,
        min_bins=3,
        min_reps=3,
        minread_threshold=0,
    )

    row = result.iloc[0]
    assert row["variant_seq"] == "variant1"
    assert np.isnan(row["Rep1.score"])
    assert np.isnan(row["Rep2.score"])
    assert np.isnan(row["Rep3.score"])
    assert np.isnan(row["score"])
