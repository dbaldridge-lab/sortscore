import pandas as pd
import pytest

from sortscore.utils.experiment_setup import load_experiment_setup


def test_load_experiment_setup_accepts_path_column(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text(
        "Replicate,Bin,Path,MFI\n"
        "1,1,counts1.tsv,123.4\n"
        "1,2,counts2.tsv,234.5\n"
    )

    df, cols = load_experiment_setup(str(setup_csv))
    assert isinstance(df, pd.DataFrame)
    assert cols.count_file == "Path"
    assert cols.replicate == "Replicate"
    assert cols.bin == "Bin"
    assert cols.mfi == "MFI"


def test_load_experiment_setup_accepts_read_counts_csv_column(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text(
        "replicate,bin,Read Counts (CSV),mfi\n"
        "1,1,counts1.tsv,123.4\n"
    )

    _, cols = load_experiment_setup(str(setup_csv))
    assert cols.count_file == "Read Counts (CSV)"
    assert cols.replicate == "replicate"
    assert cols.bin == "bin"
    assert cols.mfi == "mfi"


def test_load_experiment_setup_requires_count_file_column(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text("Replicate,Bin,MFI\n1,1,123.4\n")

    with pytest.raises(ValueError, match=r"count file path column"):
        load_experiment_setup(str(setup_csv))


def test_load_experiment_setup_validates_numeric_fields(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text(
        "Replicate,Bin,Path,MFI\n"
        "not_an_int,1,counts.tsv,123.4\n"
    )

    with pytest.raises(ValueError, match=r"Replicate column"):
        load_experiment_setup(str(setup_csv))


def test_load_experiment_setup_requires_tile_in_batch_mode(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text(
        "Replicate,Bin,Path,MFI\n"
        "1,1,counts.tsv,123.4\n"
    )

    with pytest.raises(ValueError, match=r"Tile"):
        load_experiment_setup(str(setup_csv), require_tile=True)


def test_load_experiment_setup_validates_tile_as_integer(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text(
        "Tile,Replicate,Bin,Path,MFI\n"
        "not_int,1,1,counts.tsv,123.4\n"
    )

    with pytest.raises(ValueError, match=r"Tile column"):
        load_experiment_setup(str(setup_csv), require_tile=True)


def test_load_experiment_setup_accepts_integer_tile_in_batch_mode(tmp_path):
    setup_csv = tmp_path / "experiment_setup.csv"
    setup_csv.write_text(
        "Tile,Replicate,Bin,Path,MFI\n"
        "1,1,1,counts.tsv,123.4\n"
    )

    _, cols = load_experiment_setup(str(setup_csv), require_tile=True)
    assert cols.tile == "Tile"
