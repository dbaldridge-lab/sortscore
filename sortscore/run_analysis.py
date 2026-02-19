"""
Main entry point for Sort-seq variant scoring.

This script loads the experiment configuration and orchestrates the analysis workflow.

Usage:
sortscore score -n EXPERIMENT -e path/to/experiment_setup.csv -c path/to/config.json
"""
import logging
import sys
import os
import pandas as pd
from sortscore.utils.load_experiment import ExperimentConfig
from sortscore.utils.file_utils import ensure_output_subdirs
from sortscore.utils.analysis_logger import AnalysisLogger, generate_date_suffix
from sortscore.analysis.workflows import run_variant_analysis_workflow
from sortscore.visualization.heatmap_workflow import generate_heatmap_visualizations
from sortscore.utils.console_utils import parse_analysis_args, build_merged_analysis_config


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    args = parse_analysis_args()
    
    if not args.experiment_name:
        logging.error("Missing -n/--experiment-name (required).")
        sys.exit(1)
    if not args.experiment_setup_file:
        logging.error("Missing -e/--experiment-setup-file (required).")
        sys.exit(1)

    try:
        merged, config_dir = build_merged_analysis_config(args)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(1)

    # Load experiment config using merged parameters
    try:
        experiment = ExperimentConfig.from_dict(merged, config_file_dir=config_dir)
    except Exception as e:
        logging.error(f"Failed to load experiment configuration: {e}")
        sys.exit(1)

    output_dir = experiment.output_dir or '.'
    # Ensure output directory and subdirectories exist
    ensure_output_subdirs(output_dir)

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

    # Ensure fig_format is valid
    fig_format = getattr(args, 'fig_format', None)
    if not fig_format:
        fig_format = getattr(experiment, 'fig_format', None)
    if not fig_format:
        fig_format = 'png'

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
        # Plot dispatch is controlled by mutagenesis_type.
        if experiment.mutagenesis_type == 'aa':
            if aa_scores_df is None:
                raise ValueError("mutagenesis_type='aa' requires AA scores for plotting, but none were generated.")
            generate_heatmap_visualizations(
                scores_df=aa_scores_df,
                experiment=experiment,
                output_dir=output_dir,
                output_suffix=output_suffix,
                fig_format=fig_format,
                export_positional_averages=args.pos_color,
            )
        elif experiment.mutagenesis_type == 'codon':
            if dna_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=dna_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
            if aa_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=aa_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
        else:
            # Preserve existing fallback behavior for other types (e.g., snv).
            if dna_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=dna_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
                    export_positional_averages=args.pos_color,
                )
            if aa_scores_df is not None:
                generate_heatmap_visualizations(
                    scores_df=aa_scores_df,
                    experiment=experiment,
                    output_dir=output_dir,
                    output_suffix=output_suffix,
                    fig_format=fig_format,
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
