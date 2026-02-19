"""
Main entry point for Sort-seq batch variant analysis.

This script loads batch configuration, orchestrates normalization across multiple experiments,
and generates combined outputs with proper cleanup.

Usage:
    sortscore norm --config path/to/batch_config.json
    python -m sortscore.run_batch_analysis --config path/to/batch_config.json
"""

import argparse
import json
import logging
import sys
import os
from typing import List, Optional

from sortscore.utils.load_experiment import ExperimentConfig
from sortscore.analysis.batch_config import BatchConfig
from sortscore.analysis.batch_normalization import run_batch_analysis, save_batch_results, generate_batch_visualizations


def main():
    """Main entry point for batch analysis."""
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    parser = argparse.ArgumentParser(description="Run Sort-seq batch variant analysis.")
    parser.add_argument('-c', '--config', type=str, required=True, 
                       help='Path to batch configuration JSON file')
    parser.add_argument('-s', '--suffix', type=str, 
                       help='Custom suffix for output files (default: auto-generated)')
    args = parser.parse_args()

    # Load batch configuration
    try:
        batch_config = BatchConfig.from_json(args.config)
        batch_config.validate_config()
        logging.info(f"Loaded batch config with {len(batch_config.experiment_configs)} experiments")
    except Exception as e:
        logging.error(f"Failed to load batch config: {e}")
        sys.exit(1)
    
    # Run batch analysis
    try:
        batch_config_dict = batch_config.get_batch_config_dict()
        results = run_batch_analysis(batch_config_dict)
        logging.info(f"Batch analysis complete using {results['method']} normalization")
    except Exception as e:
        logging.error(f"Failed to run batch analysis: {e}")
        sys.exit(1)
    
    # Save results
    try:
        save_batch_results(results, results['output_dir'], args.suffix)
        logging.info(f"Batch results saved to {results['output_dir']}")
    except Exception as e:
        logging.error(f"Failed to save batch results: {e}")
        sys.exit(1)
    
    # Generate visualizations
    try:
        generate_batch_visualizations(
            results, 
            batch_config_dict, 
            results['output_dir'], 
            args.suffix
        )
        logging.info("Generated batch visualizations")
    except Exception as e:
        logging.error(f"Failed to generate visualizations: {e}")
        # Don't exit on visualization failure, just warn
    
    # Cleanup individual files if requested
    if batch_config.cleanup_individual_files:
        try:
            cleanup_individual_experiment_files(batch_config.experiment_configs)
            logging.info("Cleaned up individual experiment files")
        except Exception as e:
            logging.warning(f"Failed to cleanup individual files: {e}")
    
    print(f"Batch analysis complete! Combined results saved to {results['output_dir']}")


def cleanup_individual_experiment_files(experiment_configs: List[str]) -> None:
    """
    Clean up individual experiment output files after batch processing.
    
    Parameters
    ----------
    experiment_configs : List[str]
        List of paths to individual experiment config files
    """
    import shutil
    
    for config_path in experiment_configs:
        try:
            # Load experiment config to get output directory
            experiment = ExperimentConfig.from_json(config_path)
            if experiment.output_dir and os.path.exists(experiment.output_dir):
                # Remove scores and figures directories
                scores_dir = os.path.join(experiment.output_dir, 'scores')
                figures_dir = os.path.join(experiment.output_dir, 'figures')
                
                if os.path.exists(scores_dir):
                    shutil.rmtree(scores_dir)
                    logging.debug(f"Removed {scores_dir}")
                
                if os.path.exists(figures_dir):
                    shutil.rmtree(figures_dir)
                    logging.debug(f"Removed {figures_dir}")
        except Exception as e:
            logging.warning(f"Failed to cleanup files for {config_path}: {e}")


if __name__ == "__main__":
    main()
