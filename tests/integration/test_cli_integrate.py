import subprocess
import sys
from pathlib import Path

import pandas as pd


def test_sortscore_cli_integrate_lilace_writes_input_file(tmp_path, isolated_runtime_env):
    input_path = tmp_path / "batch_scores.csv"
    output_path = tmp_path / "Lilace_input_counts.csv"

    pd.DataFrame(
        [
            {
                "batch": "tile4",
                "aa_seq_diff": "p.Gly1Ala",
                "annotate_aa": "G1A",
                "count.r1b1": 10,
                "count.r1b2": 20,
                "count.r1b3": 30,
                "count.r1b4": 40,
                "avgscore": 1.5,
            },
        ]
    ).to_csv(input_path, index=False)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "sortscore",
            "integrate",
            "lilace",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--batch",
            "tile4",
        ],
        capture_output=True,
        text=True,
        env=isolated_runtime_env,
    )

    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    assert output_path.exists()

    exported = pd.read_csv(output_path)
    assert list(exported.columns) == [
        "variant_id",
        "mutation_type",
        "position",
        "replicate",
        "c_0",
        "c_1",
        "c_2",
        "c_3",
        "sortscore_score",
        "sortscore_score_column",
    ]
