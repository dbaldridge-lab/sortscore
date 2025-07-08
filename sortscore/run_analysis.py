"""
Main entry point for Sort-seq variant analysis.

This script loads the experiment configuration, ensures output directories exist, and orchestrates the analysis workflow.

Usage:
    python -m sortscore.run_analysis --config path/to/your_config.json
"""
import argparse
import json
import logging
import sys
import pandas as pd
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.utils import ensure_output_subdirs


def main():
    parser = argparse.ArgumentParser(description="Run Sort-seq variant analysis.")
    parser.add_argument('-c', '--config', type=str, required=True, help='Path to experiment config JSON file')
    args = parser.parse_args()

    # Load experiment config using dataclass
    try:
        experiment = ExperimentConfig.from_json(args.config)
    except Exception as e:
        logging.error(f"Failed to load config: {e}")
        sys.exit(1)

    output_dir = experiment.output_dir or '.'
    # Ensure output subdirectories exist
    try:
        ensure_output_subdirs(output_dir)
    except PermissionError as e:
        logging.error(f"Permission denied: cannot create output directories in '{output_dir}'. {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Failed to create output directories in '{output_dir}': {e}")
        sys.exit(1)
    logging.info(f"Output directories ensured in {output_dir}")

    # Load experiment setup CSV (flexible for replicates and bins)
    try:
        setup_df = pd.read_csv(experiment.experiment_setup_file)
    except Exception as e:
        logging.error(f"Failed to load experiment setup file: {e}")
        sys.exit(1)

    # Group by replicate and bin, import counts flexibly
    counts = {}
    for _, row in setup_df.iterrows():
        rep = int(row['replicate'])
        bin_ = int(row['bin'])
        count_file = row['count_file']
        try:
            df = pd.read_csv(count_file, sep=None, engine='python')
        except Exception as e:
            logging.error(f"Failed to load count file {count_file}: {e}")
            continue
        counts.setdefault(rep, {})[bin_] = df
    logging.info(f"Loaded counts for {len(counts)} replicates.")
    # counts[rep][bin] gives the DataFrame for each replicate/bin

    print("Sort-seq analysis setup complete. Counts loaded for all replicates and bins.")

if __name__ == "__main__":
    main()
