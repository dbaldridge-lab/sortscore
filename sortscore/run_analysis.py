"""
Main entry point for Sort-seq variant scoring.

This script loads the experiment configuration and orchestrates the analysis workflow.

Usage:
    sortscore --config path/to/your_config.json
    python -m sortscore --config path/to/your_config.json
    python -m sortscore --config path/to/your_config.json --suffix custom_name
"""
import logging
import sys
import os
import pandas as pd
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.utils import ensure_output_subdirs
from sortscore.analysis.batch_workflow import run_batch_mode
from sortscore.analysis.analysis_logger import AnalysisLogger, generate_date_suffix
from sortscore.analysis.workflows import run_variant_analysis_workflow
from sortscore.visualization.heatmap_workflow import generate_heatmap_visualizations
from sortscore.console_utils import create_analysis_parser


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    parser = create_analysis_parser()
    args = parser.parse_args()
    
    # Handle batch processing mode for tiled experiments (e.g. multiple oligo pools)
    if args.batch:
        run_batch_mode(args.config, args.suffix)
        return

    # Load experiment config using dataclass
    try:
        experiment = ExperimentConfig.from_json(args.config)
    except Exception as e:
        logging.error(f"Failed to load config: {e}")
        sys.exit(1)

    output_dir = experiment.output_dir or '.'
    # Ensure output subdirectories exist
    ensure_output_subdirs(output_dir)

    # Load sequencing counts for all samples
    try:
        experiment.load_counts()
    except Exception as e:
        logging.error(f"Failed to load counts: {e}")
        sys.exit(1)
    logging.info(f"Loaded counts for {len(experiment.counts)} replicates.")
    # experiment.counts[rep][bin] returns a DataFrame for a sample

    print("Counts loaded.")
    
    # Generate string that will be used to name output files
    if args.suffix:
        output_suffix = args.suffix
        logging.info(f"Using custom output suffix: {output_suffix}")
    else:
        output_suffix = generate_date_suffix()
        logging.info(f"Using date-based output suffix: {output_suffix}")

    # Initialize analysis logger with automatic parameter extraction
    analysis_logger = AnalysisLogger(experiment, args, output_suffix, output_dir)
    
    # Run variant analysis workflow
    try:
        dna_scores_file, aa_scores_file = run_variant_analysis_workflow(
            experiment, output_dir, output_suffix, analysis_logger
        )
        
        # Load scores for visualization.
        # For DNA inputs, codon/snv analyses may also export AA-aggregated scores.
        dna_scores_df = pd.read_csv(dna_scores_file) if dna_scores_file and os.path.exists(dna_scores_file) else None
        aa_scores_df = pd.read_csv(aa_scores_file) if aa_scores_file and os.path.exists(aa_scores_file) else None
        if dna_scores_df is not None:
            logging.info(f"Loaded DNA scores for visualization: {len(dna_scores_df)} variants")
        if aa_scores_df is not None:
            logging.info(f"Loaded AA scores for visualization: {len(aa_scores_df)} variants")

        if dna_scores_df is None and aa_scores_df is None:
            raise ValueError("No scores files were generated")
            
    except Exception as e:
        analysis_logger.add_error(f"Failed to run variant analysis workflow: {e}")
        logging.error(f"Failed to run variant analysis workflow: {e}")
        sys.exit(1)
    
    # Generate visualizations
    try:
        # Generate plots for any score tables that exist.
        if dna_scores_df is not None:
            dna_plot_df = dna_scores_df.drop(columns=['aa_seq_diff', 'annotate_aa'], errors='ignore')
            generate_heatmap_visualizations(
                scores_df=dna_plot_df,
                experiment=experiment,
                output_dir=output_dir,
                output_suffix=output_suffix,
                fig_format=args.fig_format,
                export_positional_averages=args.pos_color,
            )
        if aa_scores_df is not None:
            generate_heatmap_visualizations(
                scores_df=aa_scores_df,
                experiment=experiment,
                output_dir=output_dir,
                output_suffix=output_suffix,
                fig_format=args.fig_format,
                export_positional_averages=args.pos_color,
            )
    except Exception as e:
        logging.error(f"Failed to generate visualizations: {e}")
        sys.exit(1)

    log_file = analysis_logger.finish()
    
    print(f"Analysis complete! Results saved to {output_dir}")
    print(f"Analysis log saved to {log_file}")

if __name__ == "__main__":
    main()
