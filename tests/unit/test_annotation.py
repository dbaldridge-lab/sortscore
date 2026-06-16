import pandas as pd
import pytest

from sortscore.analysis.annotation import annotate_scores_dataframe, classify_aa_variant, classify_dna_variant
from sortscore.utils.load_experiment import ExperimentConfig


def _make_preannotated_aa_config(variant_seq: list[str]) -> ExperimentConfig:
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
                    "variant_seq": variant_seq,
                    "count": [1] * len(variant_seq),
                }
            )
        }
    }
    return config


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
    assert annotated.loc[1, "annotate_dna"] == "snv"
    assert annotated.loc[1, "annotate_aa"] == "missense_aa"

    assert annotated.loc[2, "aa_seq_diff"] == "A.2.="
    assert annotated.loc[2, "dna_seq_diff"] != ""
    assert annotated.loc[2, "annotate_dna"] == "synonymous"
    assert annotated.loc[2, "annotate_aa"] == "synonymous"


def test_annotate_scores_dataframe_marks_multiple_synonymous_codons_positionally():
    wt_dna = "ATGGCCAAA"  # MAK
    synonymous = "ATGGCTAAG"  # MAK with codon changes at A2 and K3

    scores_df = pd.DataFrame({"variant_seq": [synonymous]})

    annotated = annotate_scores_dataframe(scores_df, wt_dna, mutagenesis_type="codon")

    assert annotated.loc[0, "aa_seq_diff"] == "A.2.=, K.3.="
    assert annotated.loc[0, "annotate_dna"] == "multiple_aa"
    assert annotated.loc[0, "annotate_aa"] == "multiple_aa"


def test_preannotated_aa_parse_errors_raise_value_error():
    config = _make_preannotated_aa_config(["M1M", "not-a-variant"])

    with pytest.raises(ValueError, match="Invalid pre-annotated amino-acid variant"):
        config._convert_aa_changes_to_annotations()


def test_preannotated_aa_no_change_variants_use_positional_no_change_strings():
    config = _make_preannotated_aa_config(["M1M", "A2S"])
    original_df = config.counts[1][1]

    config._convert_aa_changes_to_annotations()
    df = config.counts[1][1]

    assert df is not original_df
    assert df.loc[0, "aa_seq_diff"] == "M.1.="
    assert df.loc[0, "annotate_aa"] == "synonymous"
    assert df.loc[1, "aa_seq_diff"] == "A.2.S"
    assert df.loc[1, "annotate_aa"] == "missense_aa"
def test_classify_aa_variant_returns_multiple_for_multi_change_strings():
    assert classify_aa_variant("A.2.S, K.3.R") == "multiple_aa"


def test_annotate_scores_dataframe_marks_multiple_aa_changes_as_multiple():
    wt_dna = "ATGGCCAAA"  # MAK
    multi_missense = "ATGTCCAGA"  # MSR

    annotated = annotate_scores_dataframe(
        pd.DataFrame({"variant_seq": [multi_missense]}),
        wt_dna,
        mutagenesis_type="codon",
    )

    assert annotated.loc[0, "annotate_dna"] == "multiple_aa"
    assert annotated.loc[0, "annotate_aa"] == "multiple_aa"


def test_multiple_dna_changes_in_one_codon_do_not_force_multiple_classification():
    wt_dna = "ATGGCC"  # MA
    same_codon_multi_base_change = "ATGAGC"  # MS; GCC -> AGC changes two DNA bases but one AA position

    annotated = annotate_scores_dataframe(
        pd.DataFrame({"variant_seq": [same_codon_multi_base_change]}),
        wt_dna,
        mutagenesis_type="codon",
    )

    assert annotated.loc[0, "dna_seq_diff"] == "G.4.A, C.5.G"
    assert annotated.loc[0, "aa_seq_diff"] == "A.2.S"
    assert annotated.loc[0, "annotate_dna"] == "dinucleotide"
    assert annotated.loc[0, "annotate_aa"] == "missense_aa"


def test_classification_raises_for_unsplittable_diff_strings():
    with pytest.raises(ValueError, match="Could not parse amino acid changes"):
        classify_aa_variant("not-a-diff")

    with pytest.raises(ValueError, match="Could not parse DNA changes"):
        classify_dna_variant("bad-dna", "A.2.S")

    with pytest.raises(ValueError, match="cannot be combined with other DNA changes"):
        classify_dna_variant("=, A.3.G", "A.2.S")


def test_classify_dna_variant_uses_specific_change_counts():
    assert classify_dna_variant("A.3.G", "K.1.=",) == "synonymous"
    assert classify_dna_variant("A.3.G", "A.2.S") == "snv"
    assert classify_dna_variant("G.4.A, C.5.G", "A.2.S") == "dinucleotide"
    assert classify_dna_variant("G.4.A, C.5.G, C.6.T", "A.2.S") == "trinucleotide"
