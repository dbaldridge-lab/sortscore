from sortscore.utils.sequence_parsing import compare_to_reference


def test_compare_to_reference_returns_explicit_no_diff_marker():
    assert compare_to_reference("MA", "MA", no_change_marker="=") == "="
