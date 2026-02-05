import pandas as pd

from sortscore.analysis.annotation import annotate_scores_dataframe


def test_annotate_scores_dataframe_supports_dna_variant_type():
    wt_dna = "ATGGCC"  # MA
    missense = "ATGTCC"  # MS (A2S)
    synonymous = "ATGGCT"  # MA (A2A)

    scores_df = pd.DataFrame(
        {
            "variant_seq": [wt_dna, missense, synonymous],
            "avgscore": [1.0, 2.0, 1.5],
        }
    )

    annotated = annotate_scores_dataframe(scores_df, wt_dna, variant_type="dna")

    assert "dna_seq_diff" in annotated.columns
    assert "aa_seq_diff" in annotated.columns
    assert "annotate_dna" in annotated.columns
    assert "annotate_aa" in annotated.columns

    assert annotated.loc[0, "aa_seq_diff"] == ""
    assert annotated.loc[0, "dna_seq_diff"] == ""
    assert annotated.loc[0, "annotate_dna"] == "wt_dna"
    assert annotated.loc[0, "annotate_aa"] == "wt_dna"

    assert annotated.loc[1, "aa_seq_diff"] == "A.2.S"
    assert annotated.loc[1, "annotate_dna"] == "missense_dna"
    assert annotated.loc[1, "annotate_aa"] == "missense_aa"

    assert annotated.loc[2, "aa_seq_diff"] == ""
    assert annotated.loc[2, "dna_seq_diff"] != ""
    assert annotated.loc[2, "annotate_dna"] == "synonymous"
    assert annotated.loc[2, "annotate_aa"] == "synonymous"

