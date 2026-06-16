from sortscore.utils.sequence_parsing import (
    compare_aa_from_dna_reference,
    compare_to_reference,
)


def test_compare_to_reference_returns_explicit_no_diff_marker():
    assert compare_to_reference("MA", "MA") == "="


def test_compare_aa_from_dna_reference_marks_synonymous_positions():
    assert compare_aa_from_dna_reference("ATGGCC", "ATGGCT") == "A.2.="


def test_compare_aa_from_dna_reference_marks_multiple_synonymous_positions():
    assert (
        compare_aa_from_dna_reference("ATGGCCAAA", "ATGGCTAAG")
        == "A.2.=, K.3.="
    )
