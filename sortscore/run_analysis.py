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

    # Load all counts using the dataclass method
    try:
        experiment.load_counts()
    except Exception as e:
        logging.error(f"Failed to load counts: {e}")
        sys.exit(1)
    logging.info(f"Loaded counts for {len(experiment.counts)} replicates.")
    # experiment.counts[rep][bin] gives the DataFrame for each replicate/bin

    print("Sort-seq analysis setup complete. Counts loaded for all replicates and bins.")

if __name__ == "__main__":
    main()
