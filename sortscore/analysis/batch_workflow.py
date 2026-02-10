"""
Batch processing workflow for combining multiple Sort-seq experiments.

This module orchestrates the complete batch analysis workflow, including configuration
loading, analysis execution, result saving, visualization generation, and cleanup.
"""
import logging
import sys
from typing import Optional
from sortscore.analysis.batch_config import BatchConfig
from sortscore.analysis.batch_normalization import run_batch_analysis, save_batch_results, generate_batch_visualizations
from sortscore.utils.file_utils import cleanup_individual_experiment_files


def run_batch_mode(config_path: str, suffix: Optional[str] = None) -> None:
    """
    Run batch processing mode for combining multiple experiments.
    
    This function orchestrates the complete batch analysis workflow:
    1. Loads and validates batch configuration
    2. Runs batch analysis with cross-experiment normalization
    3. Saves combined results and statistics
    4. Generates tiled visualizations 
    5. Optionally cleans up individual experiment files
    
    Parameters
    ----------
    config_path : str
        Path to batch configuration JSON file
    suffix : Optional[str]
        Custom suffix for output files
        
    Examples
    --------
    >>> run_batch_mode('batch_config.json')
    >>> run_batch_mode('batch_config.json', suffix='final_analysis')
    """
    # Load batch configuration
    try:
        batch_config = BatchConfig.from_json(config_path)
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
        save_batch_results(results, results['output_dir'], suffix)
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
            suffix
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