from sortscore.utils.file_utils import _resolve_count_file_path


def test__resolve_count_file_path_prefers_setup_relative_path(tmp_path):
    setup_csv = tmp_path / "nested" / "experiment_setup.csv"
    setup_csv.parent.mkdir(parents=True)
    setup_csv.write_text("Replicate,Bin,Path,MFI\n1,1,counts.tsv,123.4\n")
    count_path = setup_csv.parent / "counts.tsv"
    count_path.write_text("variant_seq\tcount\nAAA\t1\n")

    resolved = _resolve_count_file_path(
        str(setup_csv),
        "counts.tsv",
        relative_path_base="setup",
    )
    assert resolved == count_path.resolve()


def test__resolve_count_file_path_uses_cwd_when_configured(tmp_path, monkeypatch):
    setup_csv = tmp_path / "nested" / "experiment_setup.csv"
    setup_csv.parent.mkdir(parents=True)
    setup_csv.write_text("Replicate,Bin,Path,MFI\n1,1,counts.tsv,123.4\n")

    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    count_path = cwd_dir / "counts.tsv"
    count_path.write_text("variant_seq\tcount\nAAA\t1\n")
    monkeypatch.chdir(cwd_dir)

    resolved = _resolve_count_file_path(
        str(setup_csv),
        "counts.tsv",
        relative_path_base="cwd",
    )
    assert resolved == count_path.resolve()


def test__resolve_count_file_path_setup_mode_does_not_fallback_to_cwd(tmp_path, monkeypatch):
    setup_csv = tmp_path / "nested" / "experiment_setup.csv"
    setup_csv.parent.mkdir(parents=True)
    setup_csv.write_text("Replicate,Bin,Path,MFI\n1,1,counts.tsv,123.4\n")

    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    count_path = cwd_dir / "counts.tsv"
    count_path.write_text("variant_seq\tcount\nAAA\t1\n")
    monkeypatch.chdir(cwd_dir)

    resolved = _resolve_count_file_path(
        str(setup_csv),
        "counts.tsv",
        relative_path_base="setup",
    )
    assert resolved != count_path.resolve()
    assert resolved == (setup_csv.parent / "counts.tsv").resolve()
