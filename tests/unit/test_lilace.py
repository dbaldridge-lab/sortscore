import pandas as pd
from pathlib import Path

from sortscore.integrations.lilace import score_table_to_lilace_csv


def test_score_table_to_lilace_csv_returns_expected_columns(tmp_path):
    input_path = tmp_path / "score_table.csv"
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
            }
        ]
    ).to_csv(input_path, index=False)

    score_table_to_lilace_csv(
        input_path,
        output_path,
        batch="tile4",
        metadata_columns={"sortscore_score": "avgscore"},
        metadata_constants={"sortscore_score_column": "avgscore"},
    )

    result = pd.read_csv(output_path)
    assert list(result.columns) == [
        "aa_seq_diff",
        "annotate_aa",
        "replicate",
        "bin1",
        "bin2",
        "bin3",
        "bin4",
        "sortscore_score",
        "sortscore_score_column",
    ]
