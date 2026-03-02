from pathlib import Path
import json

from sortscore.utils.experiment_setup import load_experiment_setup
from sortscore.utils.tile_configs import (
    get_tile_output_dirs,
    write_batch_config_json,
    build_tile_experiments,
)


def test_generate_batch_config_from_combined_experiment_setup(tmp_path):
    """
    Validate batch config generation using a repo-relative combined setup file.
    """
    setup_path = Path("demo_data/combined_experiment_setup.csv")
    assert setup_path.exists(), f"Missing expected test fixture: {setup_path}"

    setup_df, cols = load_experiment_setup(str(setup_path), require_tile=True)

    tile_output_dirs = get_tile_output_dirs(
        setup_df=setup_df,
        tile_col=cols.tile,
        base_output_dir=str(tmp_path / "batch_outputs"),
    )

    base_config = {
        "batch_normalization_method": "zscore_2pole",
        "combined_output_dir": str(tmp_path / "combined"),
        "wt_seq": "ATG",
        "min_pos": 1,
        "max_pos": 3,
        "mutagenesis_type": "codon",
        "avg_method": "rep-weighted",
    }
    experiments = build_tile_experiments(
        base_config=base_config,
        shared_setup_file=str(setup_path),
        tile_output_dirs=tile_output_dirs,
    )

    out_file = tmp_path / "batch_config.json"
    write_batch_config_json(
        base_config=base_config,
        tile_experiments=experiments,
        output_path=str(out_file),
    )

    payload = json.loads(out_file.read_text())
    assert "experiments" in payload
    assert len(payload["experiments"]) == setup_df[cols.tile].astype(int).nunique()

    tiles_in_payload = sorted(int(x["tile"]) for x in payload["experiments"])
    assert tiles_in_payload == sorted(tile_output_dirs.keys())

    for entry in payload["experiments"]:
        assert Path(entry["output_dir"]).name.startswith("tile_")
        assert entry["wt_seq"] == "ATG"
        assert entry["min_pos"] == 1
        assert entry["max_pos"] == 3
