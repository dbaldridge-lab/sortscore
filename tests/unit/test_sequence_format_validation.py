import pandas as pd
import pytest

from sortscore.utils.load_experiment import ExperimentConfig


def _counts_from_sequences(sequences):
    return {
        1: {
            1: pd.DataFrame(
                {
                    "variant_seq": sequences,
                    "count": [1] * len(sequences),
                }
            )
        }
    }


def _counts_from_sequences_custom_col(sequences, col_name):
    return {
        1: {
            1: pd.DataFrame(
                {
                    col_name: sequences,
                    "count": [1] * len(sequences),
                }
            )
        }
    }


def test_codon_mutagenesis_rejects_aa_sequences():
    config = ExperimentConfig(
        experiment_name="exp",
        experiment_setup_file="unused.csv",
        wt_seq="ATGGCC",
        mutagenesis_type="codon",
    )
    config.counts = _counts_from_sequences(["MKT", "MRT", "MNT"])

    with pytest.raises(ValueError, match="requires DNA input sequences"):
        config._validate_mutagenesis_type_against_counts()


def test_snv_mutagenesis_rejects_aa_sequences():
    config = ExperimentConfig(
        experiment_name="exp",
        experiment_setup_file="unused.csv",
        wt_seq="ATGGCC",
        mutagenesis_type="snv",
    )
    config.counts = _counts_from_sequences(["M1V", "R2C", "P3*"])

    with pytest.raises(ValueError, match="requires DNA input sequences"):
        config._validate_mutagenesis_type_against_counts()


def test_aa_mutagenesis_allows_dna_sequences():
    config = ExperimentConfig(
        experiment_name="exp",
        experiment_setup_file="unused.csv",
        wt_seq="MKT",
        mutagenesis_type="aa",
    )
    config.counts = _counts_from_sequences(["ATGGCC", "ATGACC", "ATGTCC"])

    config._validate_mutagenesis_type_against_counts()


def test_validation_supports_non_variant_seq_column_name():
    config = ExperimentConfig(
        experiment_name="exp",
        experiment_setup_file="unused.csv",
        wt_seq="ATGGCC",
        mutagenesis_type="codon",
    )
    config.counts = _counts_from_sequences_custom_col(["MKT", "MRT", "MNT"], "sequence")

    with pytest.raises(ValueError, match="requires DNA input sequences"):
        config._validate_mutagenesis_type_against_counts()
