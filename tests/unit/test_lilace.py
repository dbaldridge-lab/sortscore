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
                "score": 1.5,
            }
        ]
    ).to_csv(input_path, index=False)

    score_table_to_lilace_csv(
        input_path,
        output_path,
        batch="tile4",
        metadata_columns={"sortscore_score": "score"},
        metadata_constants={"sortscore_score_column": "score"},
    )

    result = pd.read_csv(output_path)
    assert list(result.columns) == [
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
    assert result.loc[0, "variant_id"] == "p.Gly1Ala"
    assert result.loc[0, "mutation_type"] == "G1A"
    assert result.loc[0, "position"] == 1


def test_score_table_to_lilace_csv_supports_codon_mode(tmp_path):
    input_path = tmp_path / "dna_score_table.csv"
    output_path = tmp_path / "Lilace_input_counts.csv"

    pd.DataFrame(
        [
            {
                "dna_seq_diff": "c.123A>G",
                "annotate_dna": "missense_dna",
                "count.r1b1": 10,
                "count.r1b2": 20,
                "count.r1b3": 30,
                "count.r1b4": 40,
            }
        ]
    ).to_csv(input_path, index=False)

    score_table_to_lilace_csv(
        input_path,
        output_path,
        mutagenesis_type="codon",
    )

    result = pd.read_csv(output_path)
    assert result.loc[0, "variant_id"] == "c.123A>G"
    assert result.loc[0, "mutation_type"] == "missense_dna"
    assert result.loc[0, "position"] == 123
    assert list(result.columns[:8]) == [
        "variant_id",
        "mutation_type",
        "position",
        "replicate",
        "c_0",
        "c_1",
        "c_2",
        "c_3",
    ]


def test_score_table_to_lilace_csv_filters_rows_without_position(tmp_path):
    input_path = tmp_path / "score_table.csv"
    output_path = tmp_path / "Lilace_input_counts.csv"

    pd.DataFrame(
        [
            {
                "aa_seq_diff": "=",
                "annotate_aa": "wt_dna",
                "count.r1b1": 5,
                "count.r1b2": 6,
            },
            {
                "aa_seq_diff": "A.2.S",
                "annotate_aa": "missense_aa",
                "count.r1b1": 10,
                "count.r1b2": 20,
            },
        ]
    ).to_csv(input_path, index=False)

    score_table_to_lilace_csv(input_path, output_path)

    result = pd.read_csv(output_path)
    assert result["variant_id"].tolist() == ["A.2.S"]
    assert result["position"].tolist() == [2]
