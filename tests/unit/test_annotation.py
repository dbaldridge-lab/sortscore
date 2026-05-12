import pandas as pd
import pytest

from sortscore.analysis.annotation import annotate_scores_dataframe
from sortscore.utils.load_experiment import ExperimentConfig


def test_annotate_scores_dataframe_supports_codon_mutagenesis_type():
    wt_dna = "ATGGCC"  # MA
    missense = "ATGTCC"  # MS (A2S)
    synonymous = "ATGGCT"  # MA (A2A)

    scores_df = pd.DataFrame(
        {
            "variant_seq": [wt_dna, missense, synonymous],
            "avgscore": [1.0, 2.0, 1.5],
        }
    )

    annotated = annotate_scores_dataframe(scores_df, wt_dna, mutagenesis_type="codon")

    assert "dna_seq_diff" in annotated.columns
    assert "aa_seq_diff" in annotated.columns
    assert "annotate_dna" in annotated.columns
    assert "annotate_aa" in annotated.columns

    assert annotated.loc[0, "aa_seq_diff"] == "="
    assert annotated.loc[0, "dna_seq_diff"] == "="
    assert annotated.loc[0, "annotate_dna"] == "wt_dna"
    assert annotated.loc[0, "annotate_aa"] == "wt_dna"

    assert annotated.loc[1, "aa_seq_diff"] == "A.2.S"
    assert annotated.loc[1, "annotate_dna"] == "missense_dna"
    assert annotated.loc[1, "annotate_aa"] == "missense_aa"

    assert annotated.loc[2, "aa_seq_diff"] == "="
    assert annotated.loc[2, "dna_seq_diff"] != ""
    assert annotated.loc[2, "annotate_dna"] == "synonymous"
    assert annotated.loc[2, "annotate_aa"] == "synonymous"


def test_preannotated_aa_parse_errors_raise_value_error():
    config = ExperimentConfig(
        experiment_name="exp",
        experiment_setup_file="unused.csv",
        wt_seq="ATGGCC",  # MA
        mutagenesis_type="aa",
    )
    config.counts = {
        1: {
            1: pd.DataFrame(
                {
                    "variant_seq": ["M1M", "not-a-variant"],
                    "count": [1, 1],
                }
            )
        }
    }

    with pytest.raises(ValueError, match="Invalid pre-annotated amino-acid variant"):
        config._convert_aa_changes_to_annotations()
