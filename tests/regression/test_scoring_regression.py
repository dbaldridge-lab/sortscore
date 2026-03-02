import json
from datetime import datetime
from pathlib import Path
import subprocess
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_sortscore_cli_aa_scores_regression_subset(config_dict, tmp_path, isolated_runtime_env):
    """AA score output should remain stable for an example subset of variants."""
    expected_suffix = datetime.now().strftime("%Y%m%d")
    output_root = (REPO_ROOT / "_test_outputs").resolve()
    fixture_path = REPO_ROOT / "tests" / "regression" / "fixtures" / "aa_score_regression_subset.csv"
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
        assert (
            merged[expected_col].round(3) == merged[actual_col].round(3)
        ).all(), f"Regression mismatch in '{column}':\n{merged[['aa_seq_diff', expected_col, actual_col]]}"
