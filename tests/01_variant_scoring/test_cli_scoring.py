import subprocess
import os
import sys
import tempfile
from pathlib import Path
import pytest


def test_sortscore_cli_runs_and_outputs(config_path, config_dict, cleanup_outputs):
    """Test that sortscore CLI for scoring pipeline runs and produces expected outputs."""
    data_dir = config_dict["output_dir"]
    
    # Verify that the console script was installed from setup.py
    sortscore_exe = Path(sys.executable).with_name("sortscore")
    assert sortscore_exe.exists(), f"Expected console script to exist at {sortscore_exe}"
    
    # Run the sortscore CLI with the provided config file
    with tempfile.TemporaryDirectory(prefix="sortscore_test_") as tmpdir:
        # Set environment variables to isolate matplotlib and cache directories
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        env["MPLCONFIGDIR"] = tmpdir
        env["XDG_CACHE_HOME"] = tmpdir
        env["HOME"] = tmpdir
        result = subprocess.run(
            [
                str(sortscore_exe),
                "-n",
                config_dict["experiment_name"],
                "-e",
                config_dict["experiment_setup_file"],
                "-c",
                config_path,
            ],
            capture_output=True,
            text=True,
            env=env,
        )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"

    output_root = Path(data_dir)
    assert output_root.exists(), f"Expected output directory to exist: {data_dir}"

    # Collect all output files
    files = [p for p in output_root.rglob("*") if p.is_file()]
    filenames = [p.name for p in files]

    # Define expected output files regardless of analysis type
    expected = {
        "stats_json": any(name.endswith(".json") and "_stats_" in name for name in filenames),
        "log_json": any(name.endswith(".log.json") for name in filenames),
    }

    # Add expected outputs based on variant type
    expected["aa_scores_csv"] = any(name.endswith(".csv") and "_aa_scores_" in name for name in filenames)
    expected["aa_heatmap"] = any("_aa_heatmap_" in name for name in filenames)
    expected["dna_scores_csv"] = any(name.endswith(".csv") and "_dna_scores_" in name for name in filenames)

    missing = [k for k, v in expected.items() if not v]
    if missing:
        formatted = "\n".join(sorted(str(p.relative_to(output_root)) for p in files))
        raise AssertionError(f"Missing expected outputs {missing} in {data_dir}\nFound:\n{formatted}")
