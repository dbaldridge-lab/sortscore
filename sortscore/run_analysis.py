"""
Main entry point for Sort-seq variant analysis.

This script loads the experiment configuration, ensures output directories exist, and orchestrates the analysis workflow.

Usage:
    python -m sortscore.run_analysis --config path/to/your_config.json
    python -m sortscore.run_analysis --config path/to/your_config.json --suffix custom_name
    python -m sortscore.run_analysis --batch --config path/to/batch_config.json
"""
import logging
import sys
import os
import pandas as pd
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.utils import ensure_output_subdirs
from sortscore.analysis.batch_workflow import run_batch_mode
from sortscore.analysis.score import calculate_full_activity_scores
from sortscore.analysis.data_processing import aggregate_aa_data
from sortscore.analysis.annotation import annotate_scores_dataframe
from sortscore.analysis.analysis_logger import AnalysisLogger, generate_date_suffix
from sortscore.analysis.statistics import calculate_replicate_statistics, round_score_columns, get_replicate_score_columns
from sortscore.analysis.summary_stats import calculate_summary_stats, save_summary_stats
from sortscore.analysis.aa_scores import process_and_save_aa_scores
from sortscore.sequence_parsing import translate_dna
from sortscore.visualization.heatmaps import plot_heatmap
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
    
    # Calculate activity scores using dataclass parameters
    try:
        logging.info("Calculating activity scores...")
        
        # Get merged counts DataFrame
        merged_df = experiment.get_merged_counts()
        
        scores_df = calculate_full_activity_scores(
            counts=experiment.counts,
            median_gfp=experiment.median_gfp,
            min_bins=experiment.bins_required,
            min_reps=experiment.reps_required,
            minread_threshold=experiment.minread_threshold,
            avg_method=experiment.avg_method,
            total_reads=experiment.total_reads,
            cell_prop=experiment.cell_prop,
            merged_df=merged_df,
            max_cv=experiment.max_cv
        )
        logging.info(f"Calculated scores for {len(scores_df)} variants.")
        analysis_logger.set_processing_stats(len(scores_df))
        
    except Exception as e:
        analysis_logger.add_error(f"Failed to calculate activity scores: {e}")
        logging.error(f"Failed to calculate activity scores: {e}")
        sys.exit(1)
    
    # Annotate sequences
    try:
        # Annotate the scores DataFrame with sequence information
        scores_df = annotate_scores_dataframe(scores_df, experiment.wt_seq, experiment.variant_type, experiment.transcript_id)
    except Exception as e:
        analysis_logger.add_warning(f"Failed to add annotations: {e}")
        logging.warning(f"Failed to add annotations: {e}")

    # Prepare scores DataFrame for saving
    scores_dir = os.path.join(output_dir, 'scores')
    scores_df_rounded = scores_df.copy()
    rep_score_columns = get_replicate_score_columns(scores_df_rounded)
    scores_df_rounded = calculate_replicate_statistics(scores_df_rounded, rep_score_columns)
    scores_df_rounded = round_score_columns(scores_df_rounded)
    
    # Save DNA scores if applicable
    if experiment.variant_type == 'dna':
        dna_scores_file = os.path.join(scores_dir, f"{experiment.experiment_name}_dna_scores_{output_suffix}.csv")
        scores_df_rounded.to_csv(dna_scores_file, index=False)
        logging.info(f"Saved DNA scores to {dna_scores_file}")

        analysis_logger.log_output_file(
            'dna_scores', 
            f"{experiment.experiment_name}_dna_scores_{output_suffix}.csv",
            dna_scores_file,
            variant_count=len(scores_df_rounded)
        )

    process_and_save_aa_scores(scores_df, experiment, scores_dir, output_suffix, analysis_logger)
    
    # Generate visualizations
    try:
        logging.info("Generating visualizations...")
        
        # Convert avg_method to column name format (replace hyphens with underscores)
        if experiment.avg_method == 'simple-avg':
            score_col = 'avgscore'
        else:
            score_col_suffix = experiment.avg_method.replace('-', '_')
            score_col = f'avgscore_{score_col_suffix}'
        
        figures_dir = os.path.join(output_dir, 'figures')
        
        # AA heatmap
        if 'aa_seq_diff' in scores_df.columns:
            # Use AA data directly for AA-only experiments, or aggregate for DNA experiments
            if experiment.variant_type == 'aa':
                # For AA-only data, use the data directly
                aa_data = scores_df[['aa_seq_diff', 'annotate_aa', score_col]].copy()
            else:
                # For DNA data, aggregate to AA level
                aa_data = aggregate_aa_data(scores_df, score_col)
            
            # Get WT score from the aggregated AA data for AA heatmap
            # For AA substitution heatmaps, average synonymous WT variants
            wt_score = pd.NA  # default
            wt_label = "WT Avg"
            if 'annotate_aa' in aa_data.columns:
                wt_subset = aa_data[aa_data['annotate_aa'] == 'synonymous']
                if len(wt_subset) > 0 and score_col in wt_subset.columns:
                    score_val = wt_subset[score_col].mean()
                    wt_score = float(score_val) if pd.notna(score_val) else pd.NA
                    logging.info(f"Found WT score from AA data (averaged from {len(wt_subset)} synonymous variants): {wt_score}")
            
            # Create a config object for AA heatmap with proper fields
            from types import SimpleNamespace
            aa_config = SimpleNamespace(
                experiment_name=experiment.experiment_name,
                num_aa=experiment.num_aa,
                num_positions=experiment.num_positions,
                min_pos=experiment.min_pos,
                wt_seq=translate_dna(experiment.wt_seq),
                variant_type='aa',
                mutagenesis_variants=experiment.mutagenesis_variants,
                position_type=experiment.position_type
            )
            
            # Set up colorbar with WT score if available
            tick_values = None
            tick_labels = None
            if not pd.isna(wt_score):
                # Get min/max values from data for colorbar range
                data_min = aa_data[score_col].min()
                data_max = aa_data[score_col].max()
                # Add WT score to tick marks if min/max are not NaN
                if pd.notna(data_min) and pd.notna(data_max):
                    tick_values = [data_min, wt_score, data_max]
                    tick_labels = [f'{data_min:.0f}', f'{wt_label}', f'{data_max:.0f}']
            
            aa_heatmap_file = os.path.join(figures_dir, f"aa_heatmap_{experiment.avg_method}_{output_suffix}.png")
            # Only pass wt_score if it's a valid number
            if pd.notna(wt_score):
                plot_heatmap(aa_data, score_col, aa_config, wt_score=float(wt_score),
                            tick_values=tick_values, tick_labels=tick_labels,
                            export_heatmap=True, output=figures_dir, fig_format=args.fig_format, export_matrix=True,
                            export_positional_averages=args.pos_color,
                            show_biophysical_properties=experiment.biophysical_prop,
                            suffix=output_suffix)
            else:
                plot_heatmap(aa_data, score_col, aa_config,
                            tick_values=tick_values, tick_labels=tick_labels,
                            export_heatmap=True, output=figures_dir, fig_format=args.fig_format, export_matrix=True,
                            export_positional_averages=args.pos_color,
                            show_biophysical_properties=experiment.biophysical_prop,
                            suffix=output_suffix)
            logging.info(f"Saved AA heatmap to {aa_heatmap_file}")
        
        # Codon heatmap  
        if experiment.variant_type == 'dna':
            # Get WT score for codon heatmap tick marks (use wt_dna for codon-level analysis)
            wt_score_codon = pd.NA  # default
            wt_label_codon = "WT"
            if 'annotate_dna' in scores_df.columns:
                wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
                if len(wt_subset) > 0 and score_col in wt_subset.columns:
                    # Average all wt_dna entries for codon heatmap
                    score_val = wt_subset[score_col].mean()
                    wt_score_codon = float(score_val) if pd.notna(score_val) else pd.NA
                    logging.info(f"Found WT score for codon heatmap (averaged from {len(wt_subset)} wt_dna): {wt_score_codon}")
            
            # Set up colorbar with WT score if available
            tick_values_codon = None
            tick_labels_codon = None
            if not pd.isna(wt_score_codon):
                # Get min/max values from data for colorbar range
                data_min = scores_df[score_col].min()
                data_max = scores_df[score_col].max()
                # Add WT score to tick marks if min/max are not NaN
                if pd.notna(data_min) and pd.notna(data_max):
                    tick_values_codon = [data_min, wt_score_codon, data_max]
                    tick_labels_codon = [f'{data_min:.0f}', f'{wt_label_codon}', f'{data_max:.0f}']
            
            codon_heatmap_file = os.path.join(figures_dir, f"codon_heatmap_{experiment.avg_method}_{output_suffix}.png")
            # Only pass wt_score if it's a valid number
            if pd.notna(wt_score_codon):
                plot_heatmap(scores_df, score_col, experiment, wt_score=float(wt_score_codon),
                            tick_values=tick_values_codon, tick_labels=tick_labels_codon,
                            export_heatmap=True, output=figures_dir, fig_format=args.fig_format, export_matrix=True,
                            export_positional_averages=args.pos_color,
                            show_biophysical_properties=experiment.biophysical_prop,
                            suffix=output_suffix)
            else:
                plot_heatmap(scores_df, score_col, experiment,
                            tick_values=tick_values_codon, tick_labels=tick_labels_codon,
                            export_heatmap=True, output=figures_dir, fig_format=args.fig_format, export_matrix=True,
                            export_positional_averages=args.pos_color,
                            show_biophysical_properties=experiment.biophysical_prop,
                            suffix=output_suffix)
            logging.info(f"Saved codon heatmap to {codon_heatmap_file}")
            
    except Exception as e:
        logging.error(f"Failed to generate visualizations: {e}")
        sys.exit(1)
    
    # Save summary statistics
    try:
        stats = calculate_summary_stats(scores_df, experiment)
        save_summary_stats(stats, experiment, scores_dir, output_suffix, analysis_logger)
    except Exception as e:
        analysis_logger.add_warning(f"Failed to save statistics: {e}")
        logging.warning(f"Failed to save statistics: {e}")
    
    log_file = analysis_logger.finish()
    
    print(f"Analysis complete! Results saved to {output_dir}")
    print(f"Analysis log saved to {log_file}")

if __name__ == "__main__":
    main()
