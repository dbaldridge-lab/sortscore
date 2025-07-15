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
import os
import pandas as pd
from sortscore.analysis.load_experiment import ExperimentConfig
from sortscore.analysis.utils import ensure_output_subdirs
from sortscore.analysis.score import calculate_full_activity_scores
from sortscore.analysis.data_processing import aggregate_aa_data
from sortscore.analysis.annotation import annotate_scores_dataframe, aggregate_synonymous_variants
from sortscore.sequence_parsing import translate_dna
from sortscore.visualization.plots import plot_heatmap


def main():
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
    
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
    
    # Calculate activity scores using dataclass parameters
    try:
        logging.info("Calculating activity scores...")
        scores_df = calculate_full_activity_scores(
            counts=experiment.counts,
            median_gfp=experiment.median_gfp,
            min_bins=experiment.bins_required,
            min_reps=experiment.reps_required,
            minread_threshold=experiment.minread_threshold,
            avg_method=experiment.avg_method,
            total_reads=experiment.total_reads
        )
        logging.info(f"Calculated scores for {len(scores_df)} variants.")
        logging.info(f"Score columns: {list(scores_df.columns)}")
    except Exception as e:
        logging.error(f"Failed to calculate activity scores: {e}")
        sys.exit(1)
    
    # Annotate sequences
    try:
        # Annotate the scores DataFrame with sequence information
        scores_df = annotate_scores_dataframe(scores_df, experiment.wt_seq, experiment.variant_type)
        logging.info("Added sequence annotations to scores DataFrame.")
    except Exception as e:
        logging.warning(f"Failed to add annotations: {e}")
    
    # Save results
    timestamp = pd.Timestamp.now().strftime('%Y%m%d')
    scores_dir = os.path.join(output_dir, 'scores')
    
    # Save DNA scores (full data)
    dna_scores_file = os.path.join(scores_dir, f"dna-scores_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.csv")
    scores_df.to_csv(dna_scores_file, index=False)
    logging.info(f"Saved DNA scores to {dna_scores_file}")
    
    # Save AA scores (aggregates synonymous variants)
    if 'aa_seq_diff' in scores_df.columns:
        # Filter out rows with NaN values first
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        scores_df_drop_nan = scores_df.dropna(subset=[score_col])
        
        # Aggregate synonymous variants
        columns_to_average = ['avgscore', 'avgscore_rep_weighted', 'avgscore_codon_weighted', 'Rep1.score', 'Rep2.score', 'Rep3.score']
        aa_scores = scores_df_drop_nan.groupby(['aa_seq_diff', 'annotate_aa'])[columns_to_average].mean().reset_index()
        
        aa_scores_file = os.path.join(scores_dir, f"aa-scores_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.csv")
        aa_scores.to_csv(aa_scores_file, index=False)
        logging.info(f"Saved aggregated AA scores to {aa_scores_file} ({len(aa_scores)} unique AA variants)")
    
    # Generate visualizations
    try:
        logging.info("Generating visualizations...")
        
        # Convert avg_method to column name format (replace hyphens with underscores)
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        
        figures_dir = os.path.join(output_dir, 'figures')
        
        # AA heatmap
        if 'aa_seq_diff' in scores_df.columns:
            # Aggregate data to AA level
            aa_data = aggregate_aa_data(scores_df, score_col)
            
            # Get WT score from the data
            wt_score = pd.NA  # default
            if 'annotate_aa' in scores_df.columns:
                wt_subset = scores_df[scores_df['annotate_aa'] == 'wt_dna']
                if len(wt_subset) > 0 and score_col in wt_subset.columns:
                    wt_score = float(wt_subset[score_col].iloc[0])
                    logging.info(f"Found WT score: {wt_score}")
            
            # Create a simple config object for AA heatmap
            from types import SimpleNamespace
            aa_config = SimpleNamespace(
                num_aa=experiment.num_aa,
                min_pos=experiment.min_pos,
                wt_seq=translate_dna(experiment.wt_seq),
                variant_type='aa'
            )
            
            # Set up colorbar with WT score if available
            tick_values = None
            tick_labels = None
            if not pd.isna(wt_score):
                # Get min/max values from data for colorbar range
                data_min = aa_data[score_col].min()
                data_max = aa_data[score_col].max()
                # Add WT score to tick marks
                tick_values = [data_min, wt_score, data_max]
                tick_labels = [f'{data_min:.0f}', f'WT={wt_score:.0f}', f'{data_max:.0f}']
            
            aa_heatmap_file = os.path.join(figures_dir, f"aa_heatmap_{experiment.avg_method}_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.png")
            plot_heatmap(aa_data, score_col, aa_config, wt_score=wt_score,
                        tick_values=tick_values, tick_labels=tick_labels,
                        export=True, output=aa_heatmap_file, format='png', export_matrix=True)
            logging.info(f"Saved AA heatmap to {aa_heatmap_file}")
        
        # Codon heatmap  
        if experiment.variant_type == 'dna':
            codon_heatmap_file = os.path.join(figures_dir, f"codon_heatmap_{experiment.avg_method}_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.png")
            plot_heatmap(scores_df, score_col, experiment,
                        export=True, output=codon_heatmap_file, format='png', export_matrix=True)
            logging.info(f"Saved codon heatmap to {codon_heatmap_file}")
            
    except Exception as e:
        logging.error(f"Failed to generate visualizations: {e}")
        sys.exit(1)
    
    # Save summary statistics
    try:
        stats = {}
        score_col_suffix = experiment.avg_method.replace('-', '_')
        score_col = f'avgscore_{score_col_suffix}'
        if score_col in scores_df.columns:
            stats['all_avg'] = round(float(scores_df[score_col].mean()))
            stats['all_min'] = round(float(scores_df[score_col].min()))
            stats['all_max'] = round(float(scores_df[score_col].max()))
            
            # Add annotation-based stats
            if 'annotate_dna' in scores_df.columns:
                # WT stats from DNA level
                wt_subset = scores_df[scores_df['annotate_dna'] == 'wt_dna']
                if len(wt_subset) > 0:
                    stats['wt_avg'] = round(float(wt_subset[score_col].mean()))
                    stats['wt_min'] = round(float(wt_subset[score_col].min()))
                    stats['wt_max'] = round(float(wt_subset[score_col].max()))
                
                # Synonymous (WT) stats from DNA level
                syn_subset = scores_df[scores_df['annotate_dna'] == 'synonymous']
                if len(syn_subset) > 0:
                    stats['synonymous_avg'] = round(float(syn_subset[score_col].mean()))
                    stats['synonymous_min'] = round(float(syn_subset[score_col].min()))
                    stats['synonymous_max'] = round(float(syn_subset[score_col].max()))
            
            # Missense and nonsense: use AA-aggregated scores with annotate_aa
            if 'aa_seq_diff' in scores_df.columns:
                aa_scores = aggregate_synonymous_variants(scores_df)
                if score_col in aa_scores.columns:
                    # Nonsense stats from AA level
                    nonsense_subset = aa_scores[aa_scores['annotate_aa'] == 'nonsense']
                    if len(nonsense_subset) > 0:
                        stats['nonsense_avg'] = round(float(nonsense_subset[score_col].mean()))
                        stats['nonsense_min'] = round(float(nonsense_subset[score_col].min()))
                        stats['nonsense_max'] = round(float(nonsense_subset[score_col].max()))
                    
                    # Missense stats from AA level
                    missense_subset = aa_scores[aa_scores['annotate_aa'] == 'missense_aa']
                    if len(missense_subset) > 0:
                        stats['missense_avg'] = round(float(missense_subset[score_col].mean()))
                        stats['missense_min'] = round(float(missense_subset[score_col].min()))
                        stats['missense_max'] = round(float(missense_subset[score_col].max()))
        
        stats_file = os.path.join(output_dir, f"stats_{experiment.avg_method}_{experiment.submission}_{experiment.bins_required}-bins_{int(experiment.minread_threshold)}-minreads_{timestamp}.json")
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        logging.info(f"Saved statistics to {stats_file}")
        
    except Exception as e:
        logging.warning(f"Failed to save statistics: {e}")
    
    print(f"Analysis complete! Results saved to {output_dir}")

if __name__ == "__main__":
    main()
