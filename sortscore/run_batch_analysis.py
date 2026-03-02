"""
Main entry point for Sort-seq batch variant analysis.

This script loads batch configuration, orchestrates normalization across multiple experiments,
and generates combined outputs.

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
    """Main entry point for batch analysis."""
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    parser = argparse.ArgumentParser(description="Run Sort-seq batch variant analysis.")
    parser.add_argument('-c', '--config', type=str, required=True, 
                       help='Path to batch configuration JSON file')
    parser.add_argument('-o', '--output-dir', type=str,
                       help='Override combined output directory from batch config JSON')
    parser.add_argument('--method', type=str, choices=['zscore_2pole', '2pole', 'zscore_center'],
                       help='Override normalization method from batch config JSON')
    parser.add_argument('-s', '--suffix', type=str, 
                       help='Custom suffix for output files (default: auto-generated)')
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
    
    print(f"Batch analysis complete! Combined results saved to {results['output_dir']}")


if __name__ == "__main__":
    main()
