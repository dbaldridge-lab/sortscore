"""
Main entry point for Sort-seq cross-tile normalization.

This script loads the batch configuration, orchestrates normalization across
multiple tile outputs, and generates combined results.

Usage:
    sortscore norm --config path/to/batch_config.json
    python -m sortscore.run_batch_analysis --config path/to/batch_config.json
"""

import argparse
import logging
import sys
import os

from sortscore.analysis.batch_config import BatchConfig
from sortscore.analysis.batch_normalization import run_batch_analysis, save_batch_results, generate_batch_visualizations


def main():
    """Main entry point for cross-tile normalization."""
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    parser = argparse.ArgumentParser(description="Run Sort-seq cross-tile normalization.")
    parser.add_argument('-c', '--config', type=str, required=True, 
                       help='Path to batch configuration JSON file')
    parser.add_argument('-o', '--output-dir', type=str,
                       help='Override combined output directory from batch config JSON')
    parser.add_argument('--method', type=str, choices=['zscore_2pole', 'linear_range', 'zscore_onepole'],
                       help='Override normalization method from batch config JSON')
    args = parser.parse_args()

    # Load batch configuration
    try:
        batch_config = BatchConfig.from_json(args.config)
        if args.method:
            batch_config.batch_normalization_method = args.method
        if args.output_dir:
            batch_config.combined_output_dir = os.path.abspath(os.path.expanduser(args.output_dir))
        batch_config.validate_config()
        logging.info(f"Loaded batch config with {len(batch_config.experiments)} tiles")
    except Exception as e:
        logging.error(f"Failed to load batch config: {e}")
        sys.exit(1)
    
    # Run batch analysis
    try:
        batch_config_dict = batch_config.get_batch_config_dict()
        results = run_batch_analysis(batch_config_dict)
        logging.info(f"Cross-tile normalization complete using {results['method']}")
    except Exception as e:
        logging.error(f"Failed to run batch analysis: {e}")
        sys.exit(1)
    
    # Save results
    try:
        save_batch_results(results, results['output_dir'])
        logging.info(f"Cross-tile results saved to {results['output_dir']}")
    except Exception as e:
        logging.error(f"Failed to save cross-tile results: {e}")
        sys.exit(1)
    
    # Generate visualizations
    try:
        generate_batch_visualizations(
            results, 
            batch_config_dict, 
            results['output_dir']
        )
        logging.info("Generated batch visualizations")
    except Exception as e:
        logging.error(f"Failed to generate visualizations: {e}")
        # Don't exit on visualization failure, just warn
    
    print(f"Cross-tile normalization complete! Results saved to {results['output_dir']}")


if __name__ == "__main__":
    main()
