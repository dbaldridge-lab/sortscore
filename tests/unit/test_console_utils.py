import json

from sortscore.utils.console_utils import build_merged_analysis_config, create_analysis_parser


def test_build_merged_analysis_config_includes_relative_path_base_from_cli():
    parser = create_analysis_parser()
    args = parser.parse_args(
        [
            "-n",
            "exp",
            "-e",
            "experiment_setup.csv",
            "--relative-path-base",
            "cwd",
        ]
    )

    merged, _ = build_merged_analysis_config(args)
    assert merged["relative_path_base"] == "cwd"


def test_build_merged_analysis_config_cli_overrides_config_relative_path_base(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "wt_seq": "ATGGCC",
                "relative_path_base": "setup",
            }
        )
    )

    parser = create_analysis_parser()
    args = parser.parse_args(
        [
            "-n",
            "exp",
            "-e",
            "experiment_setup.csv",
            "-c",
            str(config_path),
            "--relative-path-base",
            "cwd",
        ]
    )

    merged, _ = build_merged_analysis_config(args)
    assert merged["relative_path_base"] == "cwd"
