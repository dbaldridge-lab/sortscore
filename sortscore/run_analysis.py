"""
Main entry point for Sort-seq variant scoring.

This script loads the experiment configuration and orchestrates the analysis workflow.

Usage:
sortscore -n EXPERIMENT -e path/to/experiment_setup.csv -c path/to/config.json
"""
import logging
import sys
import os
import pandas as pd
import json
from pathlib import Path
from sortscore.utils.load_experiment import ExperimentConfig
from sortscore.utils.file_utils import ensure_output_subdirs
from sortscore.analysis.batch_workflow import run_batch_mode
from sortscore.utils.analysis_logger import AnalysisLogger, generate_date_suffix
from sortscore.analysis.workflows import run_variant_analysis_workflow
from sortscore.visualization.heatmap_workflow import generate_heatmap_visualizations
from sortscore.utils.console_utils import create_analysis_parser


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
    parser = create_analysis_parser()
    args = parser.parse_args()
    
    # Handle batch processing mode for tiled experiments (e.g. multiple oligo pools)
    if args.batch:
        if not args.config:
            logging.error("Batch mode requires -c/--config to be set to a batch config JSON file.")
            sys.exit(1)
        run_batch_mode(args.config, args.suffix)
        return

    if not args.experiment_name:
        logging.error("Missing -n/--experiment-name (required).")
        sys.exit(1)
    if not args.experiment_setup_file:
        logging.error("Missing -e/--experiment-setup-file (required).")
        sys.exit(1)

    # Merge config file (optional) with CLI overrides (highest precedence).
    merged = {}
    config_dir = None
    if args.config:
        try:
            config_path = Path(args.config).expanduser().resolve()
            config_dir = config_path.parent
            with open(config_path, "r") as f:
                merged.update(json.load(f))
        except Exception as e:
            logging.error(f"Failed to load config JSON '{args.config}': {e}")
            sys.exit(1)

    merged["experiment_name"] = args.experiment_name
    merged["experiment_setup_file"] = str(Path(args.experiment_setup_file).expanduser().resolve())
    if args.wt_seq is not None:
        merged["wt_seq"] = args.wt_seq
    if args.output_dir is not None:
        merged["output_dir"] = str(Path(args.output_dir).expanduser().resolve())
    if args.bins_required is not None:
        merged["bins_required"] = args.bins_required
    if args.reps_required is not None:
        merged["reps_required"] = args.reps_required
    if args.avg_method is not None:
        merged["avg_method"] = args.avg_method
    if args.minread_threshold is not None:
        merged["minread_threshold"] = args.minread_threshold
    if args.max_cv is not None:
        merged["max_cv"] = args.max_cv
    if args.mutagenesis_variants is not None:
        merged["mutagenesis_variants"] = [v.strip() for v in args.mutagenesis_variants.split(",") if v.strip()]
    if args.position_offset is not None:
        merged["position_offset"] = args.position_offset
    if args.biophysical_prop is not None:
        merged["biophysical_prop"] = bool(args.biophysical_prop)
    if args.position_type is not None:
        merged["position_type"] = args.position_type
    if args.min_pos is not None:
        merged["min_pos"] = args.min_pos
    if args.max_pos is not None:
        merged["max_pos"] = args.max_pos

    # Load experiment config using merged parameters
    try:
        experiment = ExperimentConfig.from_dict(merged, config_file_dir=config_dir)
    except Exception as e:
        logging.error(f"Failed to load experiment configuration: {e}")
        sys.exit(1)

    output_dir = experiment.output_dir or '.'
    # Ensure output subdirectories exist
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
