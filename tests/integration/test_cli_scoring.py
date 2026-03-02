import subprocess
import tempfile
import sys
from pathlib import Path
import json
from datetime import datetime
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_sortscore_cli_runs_and_outputs(config_dict, cleanup_outputs, isolated_runtime_env):
    """Test scoring CLI run produces fresh outputs for this invocation."""
    expected_suffix = datetime.now().strftime("%Y%m%d")
    output_root = (REPO_ROOT / "_test_outputs").resolve()
    pre_existing = set(output_root.rglob("test_experiment_*"))
    with tempfile.TemporaryDirectory(prefix="sortscore_cfg_") as cfg_tmp:
        run_cfg = dict(config_dict)
        run_cfg["output_dir"] = str(output_root)
        config_path = Path(cfg_tmp) / "config.json"
        config_path.write_text(json.dumps(run_cfg))

        # Run the sortscore CLI via module entrypoint to avoid installer coupling in tests.
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "sortscore",
                "score",
                "-n",
                "test_experiment",
                "-e",
                "demo_data/GLI2_oPool5b/experiment_setup.csv",
                "-c",
                str(config_path),
            ],
            capture_output=True,
            text=True,
            env=isolated_runtime_env,
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        assert output_root.exists(), f"Expected output directory to exist: {output_root}"
        created = set(output_root.rglob("test_experiment_*")) - pre_existing
        expected_files = [
            output_root / "scores" / f"test_experiment_dna_scores_{expected_suffix}.csv",
            output_root / "scores" / f"test_experiment_aa_scores_{expected_suffix}.csv",
            output_root / "scores" / f"test_experiment_dna_stats_{expected_suffix}.json",
            output_root / "scores" / f"test_experiment_aa_stats_{expected_suffix}.json",
            output_root / "figures" / f"test_experiment_aa_heatmap_{expected_suffix}.png",
            output_root / "figures" / f"test_experiment_aa_heatmap_matrix_{expected_suffix}.csv",
        ]
        missing = [str(p) for p in expected_files if not p.exists()]
        assert not missing, f"Missing expected fresh output files:\n" + "\n".join(missing)
        assert created, "Expected scoring run to create output artifacts"


def test_sortscore_cli_multitile_writes_tile_scores(config_dict, batch_config_dict, isolated_runtime_env):
    """Tiled scoring should emit per-tile score outputs."""
    expected_suffix = datetime.now().strftime("%Y%m%d")
    workdir = (REPO_ROOT / "_test_outputs").resolve()
    with tempfile.TemporaryDirectory(prefix="sortscore_multitile_cfg_") as cfg_tmp:
        output_root = workdir
        config_path = Path(cfg_tmp) / "config.json"
        setup_path = (REPO_ROOT / "demo_data" / "combined_experiment_setup.csv").resolve()

        tile_dirs = {
            int(entry["tile"]): Path(entry["output_dir"]).resolve()
            for entry in batch_config_dict["experiments"]
        }
        expected_tile_files = {
            tile: tile_dirs[tile] / "scores" / f"test_multitile_tile{tile}_aa_scores_{expected_suffix}.csv"
            for tile in (1, 2)
        }
        for p in expected_tile_files.values():
            if p.exists():
                p.unlink()

        run_cfg = dict(config_dict)
        run_cfg["experiments"] = batch_config_dict["experiments"]
        config_path.write_text(json.dumps(run_cfg))

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "sortscore",
                "score",
                "-n",
                "test_multitile",
                "-e",
                str(setup_path),
                "-c",
                str(config_path),
            ],
            capture_output=True,
            text=True,
            env=isolated_runtime_env,
        )

        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        for tile in (1, 2):
            tile_scores = expected_tile_files[tile]
            assert tile_scores.exists(), f"Missing tile score file: {tile_scores}"
            tile_fig = tile_dirs[tile] / "figures" / f"test_multitile_tile{tile}_aa_heatmap_{expected_suffix}.png"
            assert tile_fig.exists(), f"Missing tile figure file: {tile_fig}"

        # Distinct tile inputs should not produce identical score tables.
        tile1_scores = expected_tile_files[1]
        tile2_scores = expected_tile_files[2]
        df1 = pd.read_csv(tile1_scores)
        df2 = pd.read_csv(tile2_scores)
        assert not df1.equals(df2), "Tile 1 and Tile 2 score outputs are unexpectedly identical"
