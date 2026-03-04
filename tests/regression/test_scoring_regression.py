import json
from datetime import datetime
from pathlib import Path
import subprocess
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SINGLE_EXPERIMENT_FIXTURE = (
    REPO_ROOT
    / "tests"
    / "regression"
    / "fixtures"
    / "aa_score_regression_subset_cell_prop_normalization.csv"
)
MULTITILE_TILE2_FIXTURE = (
    REPO_ROOT
    / "tests"
    / "regression"
    / "fixtures"
    / "aa_score_regression_subset_tile2.csv"
)


def _assert_aa_regression_subset(output_path, fixture_path):
    assert output_path.exists(), f"Missing AA scores output file: {output_path}"
    expected_df = pd.read_csv(fixture_path)
    actual_df = pd.read_csv(output_path)
    actual_subset = actual_df[["aa_seq_diff", "avgscore", "avgscore_rep_weighted"]]
    merged = expected_df.merge(actual_subset, on="aa_seq_diff", how="left", suffixes=("_expected", "_actual"))

    missing_variants = merged[merged["avgscore_actual"].isna()]["aa_seq_diff"].tolist()
    assert not missing_variants, f"Missing example variants in output: {missing_variants}"

    for column in ["avgscore", "avgscore_rep_weighted"]:
        expected_col = f"{column}_expected"
        actual_col = f"{column}_actual"
        mismatch = merged[merged[expected_col].round(3) != merged[actual_col].round(3)][
            ["aa_seq_diff", expected_col, actual_col]
        ]
        if not mismatch.empty:
            raise AssertionError(
                f"Regression mismatch in '{column}' (showing up to 5 rows):\n"
                f"{mismatch.head(5).to_string(index=False)}"
            )


def test_sortscore_cli_aa_scores_regression_subset(config_dict, tmp_path, isolated_runtime_env):
    """AA score output should remain stable for an example subset of variants."""
    expected_suffix = datetime.now().strftime("%Y%m%d")
    output_root = (REPO_ROOT / "_test_outputs").resolve()
    output_path = output_root / "scores" / f"test_experiment_regression_aa_scores_{expected_suffix}.csv"
    if output_path.exists():
        output_path.unlink()

    run_cfg = dict(config_dict)
    run_cfg["output_dir"] = str(output_root)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(run_cfg))

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "sortscore",
            "score",
            "-n",
            "test_experiment_regression",
            "-e",
            "demo_data/GLI2_oPool5b/experiment_setup.csv",
            "-c",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env=isolated_runtime_env,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"

    _assert_aa_regression_subset(output_path, SINGLE_EXPERIMENT_FIXTURE)


def test_sortscore_cli_multitile_tile2_regression_subset(
    config_dict, batch_config_dict, tmp_path, isolated_runtime_env
):
    """Multitile scoring should preserve the expected AA scores on tile 2."""
    expected_suffix = datetime.now().strftime("%Y%m%d")
    setup_path = REPO_ROOT / "demo_data" / "combined_experiment_setup.csv"
    tile2_entry = next(entry for entry in batch_config_dict["experiments"] if int(entry["tile"]) == 2)
    tile2_output_dir = Path(tile2_entry["output_dir"]).resolve()
    output_path = tile2_output_dir / "scores" / f"test_multitile_regression_tile2_aa_scores_{expected_suffix}.csv"
    if output_path.exists():
        output_path.unlink()

    run_cfg = dict(config_dict)
    run_cfg["experiments"] = batch_config_dict["experiments"]
    config_path = tmp_path / "multitile_config.json"
    config_path.write_text(json.dumps(run_cfg))

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "sortscore",
            "score",
            "-n",
            "test_multitile_regression",
            "-e",
            str(setup_path),
            "-c",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env=isolated_runtime_env,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    _assert_aa_regression_subset(output_path, MULTITILE_TILE2_FIXTURE)
